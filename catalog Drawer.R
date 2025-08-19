# ============================
# Fetch the real IceCube 110-source catalog (PRL 2020)
# ============================

# Packages
req_pkgs <- c("fs","glue","cli","progress","httr2","rvest","xml2","readr","stringr","janitor","dplyr","purrr","tidyr")
to_install <- setdiff(req_pkgs, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(req_pkgs, library, character.only = TRUE))

# Output dirs
OUTDIR <- "outputs/catalog"
TMPDIR <- file.path(OUTDIR, "tmp")
fs::dir_create(OUTDIR, recurse = TRUE)
fs::dir_create(TMPDIR, recurse = TRUE)

# Helpers ---------------------------------------------------------------

is_probably_catalog <- function(x) {
  xlc <- tolower(x)
  any(str_detect(xlc, c("catalog","catalogue","source","table","supplement")))
}

# Parse RA/Dec from various common formats into degrees
# Accepts decimal degrees, or sexagesimal "hh:mm:ss" (RA) / "dd:mm:ss" (Dec).
sexagesimal_to_deg <- function(x, is_ra = FALSE){
  x <- str_trim(x)
  if (is.na(x) || x == "") return(NA_real_)
  # if numeric, assume already deg
  if (suppressWarnings(!is.na(as.numeric(x)))) return(as.numeric(x))
  # allow separators h m s / : / space
  parts <- str_split(x, "[hms:\\s]+", simplify = TRUE)
  parts <- parts[parts != ""]
  if (length(parts) < 3) return(NA_real_)
  h <- as.numeric(parts[1]); m <- as.numeric(parts[2]); s <- as.numeric(parts[3])
  if (any(is.na(c(h,m,s)))) return(NA_real_)
  val <- abs(h) + m/60 + s/3600
  if (is_ra) {
    # hours -> deg
    val <- val * 15
    return(val)
  } else {
    # handle sign from first field if present
    sign <- ifelse(str_detect(x, "^-"), -1, 1)
    return(sign * val)
  }
}

# Try reading any tabular file and standardize columns to {name, ra_deg, dec_deg}
standardize_catalog <- function(path){
  # try CSV/TSV/fwf via readr; fallback to space-delim
  df <- tryCatch({
    ext <- tolower(fs::path_ext(path))
    if (ext %in% c("csv")) readr::read_csv(path, show_col_types = FALSE)
    else if (ext %in% c("tsv","tab")) readr::read_tsv(path, show_col_types = FALSE)
    else if (ext %in% c("txt","dat","list")) readr::read_table(path, show_col_types = FALSE, comment = "#")
    else readr::read_delim(path, delim = ",", show_col_types = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(df) || !nrow(df)) return(NULL)
  
  df <- janitor::clean_names(df)
  
  # attempt to identify columns
  nm <- names(df)
  
  # Candidate name columns
  name_col <- nm[str_detect(nm, "name|src|source|object|id|designation")][1]
  
  # Candidate RA/Dec columns
  ra_col   <- nm[str_detect(nm, "\\bra\\b|right_?asc|alpha")][1]
  dec_col  <- nm[str_detect(nm, "dec|decl|delta")][1]
  
  # Some supplements split sexagesimal into 3 columns each: ra_h ra_m ra_s / dec_d dec_m dec_s
  rah <- nm[str_detect(nm, "^ra[_\\s]*h$|ra_h|rah")]
  ram <- nm[str_detect(nm, "^ra[_\\s]*m$|ra_m|ram")]
  ras <- nm[str_detect(nm, "^ra[_\\s]*s$|ra_s|ras")]
  dd  <- nm[str_detect(nm, "^dec[_\\s]*d$|dec_d|decd")]
  dm  <- nm[str_detect(nm, "^dec[_\\s]*m$|dec_m|decm")]
  ds  <- nm[str_detect(nm, "^dec[_\\s]*s$|dec_s|decs")]
  
  if (!is.na(ra_col) && !is.na(dec_col)) {
    out <- df %>%
      mutate(
        ra_deg  = purrr::map_dbl(.data[[ra_col]], ~sexagesimal_to_deg(as.character(.x), is_ra = TRUE)),
        dec_deg = purrr::map_dbl(.data[[dec_col]], ~sexagesimal_to_deg(as.character(.x), is_ra = FALSE))
      )
  } else if (all(c(length(rah), length(ram), length(ras), length(dd), length(dm), length(ds)) >= 1)) {
    # combine sexagesimal components
    out <- df %>%
      mutate(
        ra_str  = sprintf("%s:%s:%s", .data[[rah[1]]], .data[[ram[1]]], .data[[ras[1]]]),
        dec_str = sprintf("%s:%s:%s", .data[[dd[1]]],  .data[[dm[1]]],  .data[[ds[1]]]),
        ra_deg  = purrr::map_dbl(ra_str,  ~sexagesimal_to_deg(.x, is_ra = TRUE)),
        dec_deg = purrr::map_dbl(dec_str, ~sexagesimal_to_deg(.x, is_ra = FALSE))
      )
  } else {
    # Heuristics: look for columns literally named ra/dec in degrees
    if ("ra" %in% nm && "dec" %in% nm) {
      out <- df %>% mutate(ra_deg = as.numeric(ra), dec_deg = as.numeric(dec))
    } else {
      return(NULL)
    }
  }
  
  # Name fallback
  if (is.na(name_col)) {
    # look for a common catalog name field
    cand <- nm[str_detect(nm, "ngc|pks|txs|3fgl|4fgl|bll|blazar|messier|gb6")]
    name_col <- if (length(cand)) cand[1] else nm[1]
  }
  
  out %>%
    transmute(
      name = as.character(.data[[name_col]]),
      ra_deg = as.numeric(ra_deg),
      dec_deg = as.numeric(dec_deg)
    ) %>%
    filter(!is.na(ra_deg), !is.na(dec_deg), !is.na(name)) %>%
    distinct()
}

# Downloader with progress
download_to <- function(url, dest){
  cli::cli_inform(glue("Downloading: {url}"))
  resp <- httr2::request(url) |> httr2::req_user_agent("R/httr2 IceCube catalog fetcher") |> httr2::req_perform()
  if (resp$status_code >= 400) stop(glue("HTTP {resp$status_code} for {url}"))
  writeBin(httr2::resp_body_raw(resp), dest)
  return(dest)
}

# 1) Primary: APS Supplemental (PRL 124, 051103) ------------------------

aps_supp_url <- "https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.124.051103"

cli::cli_h2("Step 1: Try APS Supplemental (PRL 124, 051103)")
cat_files <- character(0)

try({
  pg <- httr2::request(aps_supp_url) |> httr2::req_perform() |> httr2::resp_body_string()
  html <- read_html(pg)
  links <- html |>
    rvest::html_elements("a") |>
    rvest::html_attr("href") |>
    unique() |>
    na.omit()
  
  # Keep only APS supplemental links (absolute)
  links <- links[str_detect(links, "^https?://")]
  
  # Narrow to likely data/tabular payloads
  links_keep <- links[str_detect(tolower(links), "\\.(csv|tsv|txt|zip|tar|gz)$") | is_probably_catalog(basename(links))]
  
  if (length(links_keep)) {
    pb <- progress_bar$new(total = length(links_keep), format = "  downloading [:bar] :current/:total eta: :eta")
    for (u in links_keep) {
      pb$tick()
      dest <- file.path(TMPDIR, basename(u))
      download_to(u, dest)
      if (str_ends(u, regex("\\.zip$", ignore_case = TRUE))) {
        unzip(dest, exdir = TMPDIR)
      } else if (str_ends(u, regex("\\.(tar|tar\\.gz|tgz|gz)$", ignore_case = TRUE))) {
        # best effort untar
        try(utils::untar(dest, exdir = TMPDIR), silent = TRUE)
      }
    }
  }
  cat_files <- fs::dir_ls(TMPDIR, recurse = TRUE, type = "file")
}, silent = TRUE)

# 2) Fallback: Harvard Dataverse dataset for the 10-yr release ----------

if (!length(cat_files)) {
  cli::cli_h2("Step 2: Try Harvard Dataverse (IceCube data release DOI)")
  # Dataverse landing & API docs:
  # https://dataverse.harvard.edu/ ; dataset DOI: 10.7910/DVN/VKL316 (from IceCube page)
  # We'll hit "native" API to list files for the persistentId, then download all and scan.
  base_api <- "https://dataverse.harvard.edu/api"
  list_api <- glue("{base_api}/datasets/:persistentId/versions/:latest/files?persistentId=doi:10.7910/DVN/VKL316")
  
  # list files JSON
  resp <- try(httr2::request(list_api) |> httr2::req_perform(), silent = TRUE)
  if (!inherits(resp, "try-error") && resp$status_code < 400) {
    js <- jsonlite::fromJSON(httr2::resp_body_string(resp))
    if (!is.null(js$data) && length(js$data)) {
      # choose likely catalog-y files
      meta <- tibble::as_tibble(js$data) %>%
        mutate(label = purrr::map_chr(dataFile, ~.x$filename)) %>%
        mutate(download_url = glue("{base_api}/access/datafile/{purrr::map_chr(dataFile, ~.x$id)}"))
      
      meta_keep <- meta %>%
        filter(is_probably_catalog(label) | str_detect(tolower(label), "\\.(csv|tsv|txt|zip|tar|gz)$"))
      
      if (nrow(meta_keep)) {
        pb <- progress_bar$new(total = nrow(meta_keep), format = "  dataverse [:bar] :current/:total eta: :eta")
        for (i in seq_len(nrow(meta_keep))) {
          pb$tick()
          dest <- file.path(TMPDIR, meta_keep$label[i])
          download_to(meta_keep$download_url[i], dest)
          if (str_ends(dest, regex("\\.zip$", ignore_case = TRUE))) unzip(dest, exdir = TMPDIR)
          if (str_ends(dest, regex("\\.(tar|tar\\.gz|tgz|gz)$", ignore_case = TRUE))) {
            try(utils::untar(dest, exdir = TMPDIR), silent = TRUE)
          }
        }
        cat_files <- fs::dir_ls(TMPDIR, recurse = TRUE, type = "file")
      }
    }
  }
}

if (!length(cat_files)) {
  stop("Could not retrieve any supplemental/catalog files from APS or Dataverse. Please check your network or try again later.")
}

# 3) Parse every candidate file and merge into a clean catalog ----------

cli::cli_h2("Step 3: Parse & standardize candidate files")
cand_tab <- tibble(file = cat_files) %>%
  filter(str_detect(tolower(file), "\\.(csv|tsv|txt|dat|list)$")) %>%
  mutate(parsed = purrr::map(file, standardize_catalog)) %>%
  filter(purrr::map_lgl(parsed, ~!is.null(.x) && nrow(.x) > 0)) %>%
  mutate(n = purrr::map_int(parsed, nrow)) %>%
  arrange(desc(n))

if (!nrow(cand_tab)) stop("Downloaded files did not contain any parseable RA/Dec tables.")

# adopt the largest table as primary (expected to be the 110-source catalog)
primary <- cand_tab$parsed[[1]] %>%
  mutate(name = str_squish(name)) %>%
  distinct(name, .keep_all = TRUE)

# sanity checks: are the four headline sources present?
must_have <- c("NGC 1068","TXS 0506+056","PKS 1424+240","GB6 J1542+6129")
present <- tolower(primary$name)
missing <- must_have[!tolower(must_have) %in% present]

if (length(missing)) {
  cli::cli_warn(glue("Warning: expected sources missing from parsed table: {toString(missing)}"))
}

# Save catalog ----------------------------------------------------------

out_csv <- file.path(OUTDIR, "source_catalog_prl2020.csv")
readr::write_csv(primary, out_csv)
cli::cli_alert_success(glue("Saved catalog with {nrow(primary)} rows -> {out_csv}"))

# Quick peek
print(head(primary, 8))

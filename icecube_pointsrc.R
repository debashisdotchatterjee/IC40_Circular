# icecube_pointsrc.R (Catalog Added)
# End-to-end IceCube point-source analysis for IC40 (and later seasons)
# Unbinned space–energy likelihood, RA-scrambled post-trial calibration,
# Bayesian HMC (single-source + hierarchical stacking),
# with professional plots/tables auto-saved.

# ---------- 0) Packages & Project Setup ----------
req <- c(
  "data.table","dplyr","ggplot2","ggpointdensity","RColorBrewer","patchwork",
  "scales","readr","stringr","lubridate","purrr","tidyr","glue","magrittr",
  "here","withr","fs","cli","progress","pbapply","grid","gridExtra","cowplot",
  "matrixStats","viridis","circular","DescTools","sf","sp","Rcpp",
  "cmdstanr","posterior","bayesplot","rstan","knitr","kableExtra","xtable"
)
miss <- setdiff(req, rownames(installed.packages()))
if(length(miss)) install.packages(miss, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Prefer cmdstanr (fast, shows progress); we silently fall back to rstan if needed.
use_cmdstan <- cmdstanr::cmdstan_version(TRUE)
if (is.null(use_cmdstan)) cli::cli_alert_warning("CmdStan not found; Bayesian step will try rstan fallback.")

options(mc.cores = parallel::detectCores())

# ---------- 1) Folders ----------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- here::here(glue("pointsrc_outputs_{stamp}"))
dir_create <- function(p) if(!fs::dir_exists(p)) fs::dir_create(p, recurse = TRUE)
dirs <- c("figures","tables","diagnostics","cache","irfs")
walk(file.path(OUTDIR, dirs), dir_create)
cli::cli_h1("Outputs -> {OUTDIR}")

# ---------- 2) Helpers (angles, time, sky) ----------
deg2rad <- function(x) x * pi/180
rad2deg <- function(x) x * 180/pi
wrap360  <- function(x) (x %% 360)
clamp    <- function(x,a,b) pmax(a, pmin(b,x))

# Great-circle separation (degrees) between (RA1,Dec1) and (RA2,Dec2) in degrees
angsep_deg <- function(ra1, dec1, ra2, dec2) {
  ra1r <- deg2rad(ra1); dec1r <- deg2rad(dec1)
  ra2r <- deg2rad(ra2); dec2r <- deg2rad(dec2)
  d <- acos(pmin(1, pmax(-1, sin(dec1r)*sin(dec2r) + cos(dec1r)*cos(dec2r)*cos(ra1r-ra2r))))
  rad2deg(d)
}

# South Pole transform (IceCube): given Az,Zen (degrees) and MJD, compute RA/Dec
# Assumes lat = -90 deg, lon = 0 deg (conventionally).
# At the pole: Dec ≈ -(90 - Zenith); HourAngle H ≈ Az (approx; needs a convention).
# Compute GMST from MJD; RA = LST - H.
gmst_from_mjd <- function(mjd) {   # Vallado approx
  JD <- mjd + 2400000.5
  T <- (JD - 2451545.0)/36525
  gmst_sec <- 67310.54841 + (876600*3600 + 8640184.812866)*T + 0.093104*T^2 - 6.2e-6*T^3
  wrap360(gmst_sec/240) # degrees
}
ra_dec_from_azzen_pole <- function(az_deg, zen_deg, mjd) {
  dec <- -(90 - zen_deg)
  H   <- wrap360(az_deg)  # azimuth measured from South? If from North, use (360 - az_deg)
  lst <- gmst_from_mjd(mjd)  # longitude ~ 0 at South Pole station by convention
  ra  <- wrap360(lst - H)
  tibble(RA = ra, Dec = dec)
}

# Small-angle spherical Gaussian kernel (degrees). Uses AngErr (>= 0.2 deg).
S_space_gauss <- function(dpsi_deg, sigma_deg) {
  s <- pmax(sigma_deg, 0.2)
  (1/(2*pi*s^2)) * exp( - (dpsi_deg^2)/(2*s^2) )
}

# Von Mises-Fisher (optional exact sphere) using kappa ~ 1/sigma^2 (sigma in radians)
S_space_vmf <- function(dpsi_deg, sigma_deg) {
  srad <- pmax(deg2rad(sigma_deg), deg2rad(0.2))
  kappa <- 1/(srad^2)
  # density ∝ exp(kappa * cos(theta)); normalized on S^2 by k/(4π sinh k)
  # For small s: exp(-theta^2/(2 s^2)) ≈ exp(kappa * (1 - theta^2/2)) up to constants
  theta <- deg2rad(dpsi_deg)
  dens <- (kappa/(4*pi*sinh(kappa))) * exp(kappa * cos(theta))
  as.numeric(dens)
}

# Power-law signal density for x = log10(E/GeV)
# If E ~ E^{-gamma}, then density of x is ∝ 10^{(1 - gamma) x}. We normalize on [xmin, xmax].
S_energy_pl_log10 <- function(x, gamma, xmin, xmax) {
  if (abs(1 - gamma) < 1e-6) {
    # gamma ~ 1 => flat in lnE => flat in x; normalize to 1/(xmax - xmin)
    return(rep(1/(xmax - xmin), length(x)))
  } else {
    cst <- ( (10^((1-gamma)*x)) * log(10) ) / ( 10^((1-gamma)*xmax) - 10^((1-gamma)*xmin) )
    return(as.numeric(cst))
  }
}

# Background energy density via KDE within declination bands
make_bkg_energy_kde <- function(df, nbands = 10) {
  # df must have columns: Dec, Elog10
  qs <- seq(0,1,length.out=nbands+1)
  dec_cuts <- as.numeric(quantile(df$Dec, probs = qs, na.rm = TRUE))
  bands <- cut(df$Dec, breaks = unique(dec_cuts), include.lowest = TRUE)
  band_levels <- levels(bands)
  models <- vector("list", length(band_levels))
  for (k in seq_along(band_levels)) {
    idx <- which(bands == band_levels[k])
    x <- df$Elog10[idx]
    models[[k]] <- density(x, n = 2048, bw = "nrd0", from = min(x), to = max(x), na.rm = TRUE)
  }
  list(cuts = unique(dec_cuts), levels = band_levels, models = models)
}
eval_bkg_energy_kde <- function(elog10, dec, kde) {
  # Find band for each dec
  cuts <- kde$cuts
  # bin index
  b <- pmin(findInterval(dec, cuts, all.inside = TRUE), length(kde$models))
  dens <- numeric(length(elog10))
  for (k in seq_along(kde$models)) {
    idx <- which(b == k)
    if (length(idx)) dens[idx] <- approx(kde$models[[k]]$x, kde$models[[k]]$y, elog10[idx], rule = 2)$y
  }
  # Normalize to be proper over observed range
  dens <- pmax(dens, .Machine$double.eps)
  dens / max(1e-12, mean(dens, na.rm = TRUE)) * mean(dens, na.rm = TRUE) # stable
}

# ---------- 3) Load IC40 events ----------
infile <- "IC40_exp.csv"
stopifnot(file.exists(infile))
raw <- suppressWarnings(fread(infile))

# Standardize columns (robust to slight header variations)
names(raw) <- str_trim(names(raw))
std_name <- function(nm) {
  nm <- str_replace_all(nm, "\\[.*?\\]", "")
  nm <- str_replace_all(nm, "[^A-Za-z0-9]+", "_")
  tolower(nm)
}
names(raw) <- std_name(names(raw))

# Guess columns
findcol <- function(cands) {
  k <- which(tolower(names(raw)) %in% tolower(cands))
  if (length(k)) names(raw)[k[1]] else NA_character_
}
col_mjd   <- findcol(c("mjd","mjd_days"))
col_elog  <- findcol(c("log10_e_gev","log10_e_gev_","log10_e_geV","elog10"))
col_ang   <- findcol(c("angerr_deg","angerr","ang_err_deg"))
col_ra    <- findcol(c("ra_deg","ra"))
col_dec   <- findcol(c("dec_deg","dec"))
col_az    <- findcol(c("azimuth_deg","azimuth"))
col_zen   <- findcol(c("zenith_deg","zenith"))

df <- raw %>%
  transmute(
    MJD   = as.numeric(get(col_mjd)),
    Elog10= as.numeric(get(col_elog)),
    AngErr= as.numeric(get(col_ang)),
    RA    = wrap360(as.numeric(get(col_ra))),
    Dec   = as.numeric(get(col_dec)),
    Az    = as.numeric(get(col_az)),
    Zen   = as.numeric(get(col_zen))
  ) %>% 
  filter(is.finite(MJD), is.finite(Elog10), is.finite(AngErr), is.finite(RA), is.finite(Dec))

cli::cli_alert_success("Loaded {nrow(df)} events with columns: MJD, Elog10, AngErr, RA, Dec, Az, Zen")

# ---------- 4) Descriptives & sanity checks ----------
desc <- list(
  MJD_min = min(df$MJD), MJD_max = max(df$MJD),
  E_med   = median(df$Elog10), E_q05 = quantile(df$Elog10,.05), E_q95=quantile(df$Elog10,.95),
  S_med   = median(df$AngErr), S_q05  = quantile(df$AngErr,.05), S_q95 = quantile(df$AngErr,.95),
  frac_s1 = mean(df$AngErr < 1, na.rm=TRUE)
)
desc_tbl <- tibble(
  Metric = names(desc),
  Value  = unlist(desc)
)
write_csv(desc_tbl, file.path(OUTDIR,"tables","ic40_descriptives.csv"))

# RA uniformity (Rayleigh test)
ra_circ <- circular::circular(deg2rad(df$RA), type="angles", units="radians", template="geographics")
ray_p <- as.numeric(circular::rayleigh.test(ra_circ)$p.value)
write_csv(tibble(test="Rayleigh", p_value=ray_p), file.path(OUTDIR,"tables","ra_uniformity_test.csv"))

# ---------- 5) Pretty plots ----------
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold"),
    strip.text = element_text(face="bold"),
    legend.position = "right"
  )

p_E <- ggplot(df, aes(Elog10)) +
  geom_histogram(aes(y=after_stat(density)), bins = 80, alpha=.8) +
  geom_density(linewidth=1) +
  labs(title="Reconstructed log10(E/GeV)", x=expression(log[10](E/GeV)), y="Density") +
  theme_paper
ggsave(file.path(OUTDIR,"figures","elog10_hist_density.png"), p_E, width=7, height=4.2, dpi=220)

p_S <- ggplot(df, aes(AngErr)) +
  geom_histogram(bins=80, alpha=.8) +
  labs(title="Per-event angular error (deg)", x="AngErr [deg]", y="Count") +
  theme_paper
ggsave(file.path(OUTDIR,"figures","angerr_hist.png"), p_S, width=7, height=4.2, dpi=220)

p_sky <- ggplot(df, aes(x=RA, y=Dec)) +
  ggpointdensity::geom_pointdensity(size=0.3) +
  scale_color_viridis_c(option="plasma") +
  coord_cartesian(xlim=c(0,360), ylim=c(-90,90), expand=FALSE) +
  labs(title="Sky map: IC40 events", x="RA [deg]", y="Dec [deg]", color="Local density") +
  theme_paper
ggsave(file.path(OUTDIR,"figures","sky_density.png"), p_sky, width=7.8, height=4.6, dpi=220)

# ---------- 6) Background energy KDE (declination bands) ----------
cli::cli_h1("Fitting background energy KDE by declination bands")
set.seed(1)
bkg_kde <- make_bkg_energy_kde(df, nbands = 10)
saveRDS(bkg_kde, file.path(OUTDIR,"cache","bkg_kde_ic40.rds"))

# ---------- 7) Source catalog ----------
# If you have your own catalog CSV, point here; else we create a small demo:
catalog_file <- file.path(OUTDIR,"cache","source_catalog_demo.csv")
demo_cat <- tribble(
  ~name,           ~ra_deg,     ~dec_deg,
  "NGC 1068",      40.6698792,  -0.01329,    # STScI JWST target list
  "TXS 0506+056",  77.43,        5.72,       # Science 2018 best-fit
  "PKS 1424+240",  216.7517,    23.8000,     # 14:27:00.4 +23:48:00
  "GB6 J1542+6129",235.73725,   61.498694    # 15:42:56.94 +61:29:55.3
)
readr::write_csv(demo_cat, catalog_file)
cat_df <- readr::read_csv(catalog_file, show_col_types = FALSE)

# Mark sources on sky plot
p_sky2 <- p_sky +
  geom_point(data = cat_df, aes(ra_deg, dec_deg), color="black", shape=21, fill="yellow", size=3, stroke=0.8) +
  geom_text(data = cat_df, aes(ra_deg, dec_deg, label = name), nudge_y = 4, size=3, fontface="bold")
ggsave(file.path(OUTDIR,"figures","sky_with_sources.png"), p_sky2, width=8, height=4.8, dpi=220)

# ---------- 8) Signal & Background likelihood building blocks ----------
# Compute per-event spatial densities for a given source
compute_S_space <- function(source_ra, source_dec, angerr_deg, evRA, evDec, kernel = c("gauss","vmf")) {
  kernel <- match.arg(kernel)
  dpsi <- angsep_deg(evRA, evDec, source_ra, source_dec)
  if (kernel == "gauss") {
    S_space_gauss(dpsi, angerr_deg)
  } else {
    S_space_vmf(dpsi, angerr_deg)
  }
}

# Background space density: RA-scramble within declination bands approximates uniform RA.
# Implement as constant in RA per dec band (empirical): B_space(Dec) ∝ 1.
eval_B_space <- function(dec) rep(1, length(dec))  # uniform over RA for time-integrated

# Energy terms:
eval_S_energy <- function(elog10, dec, gamma, kde_bkg, xmin, xmax, mode = c("fallback_pl","irf")) {
  mode <- match.arg(mode)
  if (mode == "fallback_pl") {
    S <- S_energy_pl_log10(elog10, gamma, xmin, xmax)
  } else {
    stop("IRF mode not yet wired: place IRF CSVs under irfs/ and implement convolution here.")
  }
  S
}
eval_B_energy <- function(elog10, dec, kde_bkg) {
  eval_bkg_energy_kde(elog10, dec, kde_bkg)
}

# ---------- 9) Per-source frequentist fit (MLE) ----------
negloglik_source <- function(par, S_space, elog10, dec, B_space, kde_bkg, xmin, xmax, mode_energy) {
  mu     <- pmax(par[1], 0)               # non-negative
  gamma  <- clamp(par[2], 1.3, 3.8)       # sensible range
  N      <- length(elog10)
  
  S_en   <- eval_S_energy(elog10, dec, gamma, kde_bkg, xmin, xmax, mode = mode_energy)
  B_en   <- eval_B_energy(elog10, dec, kde_bkg)
  S_i    <- S_space * S_en
  B_i    <- B_space * B_en
  
  Li <- (mu/N) * S_i + (1 - mu/N) * B_i
  -sum(log(pmax(Li, 1e-300)))
}
fit_source_mle <- function(src, df, kde_bkg, mode_energy = "fallback_pl") {
  xmin <- min(df$Elog10); xmax <- max(df$Elog10)
  Ssp  <- compute_S_space(src$ra_deg, src$dec_deg, df$AngErr, df$RA, df$Dec, kernel="gauss")
  Bsp  <- eval_B_space(df$Dec)
  
  # Initial guesses
  par0 <- c(mu = 10, gamma = 2.2)
  fit  <- nlminb(par0, negloglik_source, lower = c(0,1.2), upper = c(Inf, 4.2),
                 S_space=Ssp, elog10=df$Elog10, dec=df$Dec, B_space=Bsp,
                 kde_bkg=kde_bkg, xmin=xmin, xmax=xmax, mode_energy=mode_energy)
  # logL at MLE and at mu=0 (background-only)
  nll_hat <- fit$objective
  nll0 <- negloglik_source(c(0,2.2), S_space=Ssp, elog10=df$Elog10, dec=df$Dec, B_space=Bsp,
                           kde_bkg=kde_bkg, xmin=xmin, xmax=xmax, mode_energy=mode_energy)
  TS <- 2*(nll0 - nll_hat)
  list(
    src = src$name, mu_hat = pmax(fit$par[1],0), gamma_hat = clamp(fit$par[2],1.3,3.8),
    TS = max(TS,0), conv = (fit$convergence==0), nll_hat = nll_hat, nll0 = nll0
  )
}

# ---------- 10) Post-trial calibration by RA-scrambling ----------
# We keep (Dec, Elog10, AngErr, MJD, Az, Zen) and randomize RA in [0,360) OR 
# time-scramble -> recompute RA from (Az,Zen,MJD*), which preserves az/zen non-uniformity.
scramble_ra_uniform <- function(n) runif(n, 0, 360)
scramble_time_recompute_radec <- function(df) {
  # sample MJD* from events, recompute RA from (Az,Zen,MJD*)
  mjd_star <- sample(df$MJD, size = nrow(df), replace = TRUE)
  radec <- ra_dec_from_azzen_pole(df$Az, df$Zen, mjd_star)
  radec$RA
}

calibrate_posttrial <- function(src, df, kde_bkg, n_sims = 2000, mode_energy="fallback_pl", method=c("uniformRA","time")) {
  method <- match.arg(method)
  xmin <- min(df$Elog10); xmax <- max(df$Elog10)
  Bsp  <- eval_B_space(df$Dec)
  
  pb <- progress::progress_bar$new(
    format = glue("RA-scramble {src$name} [:bar] :percent eta: :eta"),
    total = n_sims, clear = FALSE, width = 70
  )
  
  TS_null <- pbsapply(seq_len(n_sims), function(b) {
    if (method=="uniformRA") {
      RA_star <- scramble_ra_uniform(nrow(df))
    } else {
      RA_star <- scramble_time_recompute_radec(df)
    }
    Ssp <- compute_S_space(src$ra_deg, src$dec_deg, df$AngErr, RA_star, df$Dec, kernel="gauss")
    f <- nlminb(c(0,2.2), negloglik_source, lower=c(0,1.2), upper=c(Inf,4.2),
                S_space=Ssp, elog10=df$Elog10, dec=df$Dec, B_space=Bsp, kde_bkg=kde_bkg,
                xmin=xmin, xmax=xmax, mode_energy=mode_energy)
    nll0 <- negloglik_source(c(0,2.2), S_space=Ssp, elog10=df$Elog10, dec=df$Dec, B_space=Bsp,
                             kde_bkg=kde_bkg, xmin=xmin, xmax=xmax, mode_energy=mode_energy)
    TS <- 2*(nll0 - f$objective)
    pb$tick()
    max(TS,0)
  })
  TS_null
}

# ---------- 11) Run per-source frequentist analysis ----------
cli::cli_h1("Per-source MLE and post-trial calibration")
results <- list()
for (k in seq_len(nrow(cat_df))) {
  src <- cat_df[k,]
  cli::cli_alert_info("Fitting source: {src$name}")
  fit <- fit_source_mle(src, df, bkg_kde, mode_energy="fallback_pl")
  TS_null <- calibrate_posttrial(src, df, bkg_kde, n_sims = 1000, method="uniformRA")
  
  p_pre  <- 0.5* (1 - pchisq(fit$TS, df=1))  # Chernoff
  p_post <- mean(TS_null >= fit$TS)
  
  results[[k]] <- tibble(
    source = src$name, ra = src$ra_deg, dec = src$dec_deg,
    mu_hat = fit$mu_hat, gamma_hat = fit$gamma_hat, TS = fit$TS,
    p_pre = p_pre, p_post = p_post, converged = fit$conv
  )
  
  # Save TS null plot
  dnull <- tibble(TS = TS_null)
  pnull <- ggplot(dnull, aes(TS)) + geom_histogram(bins=60, fill="grey70", color="white") +
    geom_vline(xintercept = fit$TS, color="red", linewidth=1) +
    annotate("text", x = fit$TS, y = Inf, label = glue(" TS_obs = {round(fit$TS,2)}"),
             vjust = 1.5, hjust = -0.05, color="red") +
    labs(title = glue("Null TS distribution: {src$name}"), x="TS (scrambled)", y="Count") +
    theme_paper
  ggsave(file.path(OUTDIR,"figures", glue("TS_null_{str_replace_all(src$name,'[^A-Za-z0-9]+','_')}.png")),
         pnull, width=7, height=4.4, dpi=220)
  
  # Save per-source spatial weight map (Δψ vs AngErr)
  dpsi <- angsep_deg(df$RA, df$Dec, src$ra_deg, src$dec_deg)
  psp <- ggplot(tibble(dpsi=dpsi, angerr=df$AngErr), aes(angerr, dpsi)) +
    geom_hex(bins=60) + scale_fill_viridis_c() +
    labs(title=glue("Spatial proximity vs. AngErr: {src$name}"),
         x="AngErr [deg]", y=expression(Delta*psi~"[deg]")) + theme_paper
  ggsave(file.path(OUTDIR,"figures", glue("spatial_proximity_{str_replace_all(src$name,'[^A-Za-z0-9]+','_')}.png")),
         psp, width=7, height=4.6, dpi=220)
}

res_tbl <- bind_rows(results)
readr::write_csv(res_tbl, file.path(OUTDIR,"tables","per_source_results.csv"))

# Pretty LaTeX/HTML table
kt <- res_tbl %>%
  mutate(across(c(mu_hat,gamma_hat,TS,p_pre,p_post), ~round(.,4))) %>%
  arrange(p_post) %>%
  knitr::kable("html", escape = TRUE, caption = "Per-source fits and post-trial p-values") %>%
  kableExtra::kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover","condensed"))
htmltools::save_html(kt, file.path(OUTDIR,"tables","per_source_results.html"))

# ---------- 12) Bayesian HMC (single source) ----------
# We fit a single source (e.g., NGC 1068) with parameters mu >= 0 and gamma ∈ [1.3,3.8].
# Likelihood: sum log( (mu/N) * S_space * S_energy(gamma) + (1 - mu/N)*B_space*B_energy )
stan_code_single <- "
data {
  int<lower=1> N;
  vector<lower=0>[N] Sspace;
  vector<lower=0>[N] Ben;
  vector<lower=0>[N] Sen;
}
parameters {
  real<lower=0> mu;          // expected signal count
  real<lower=1.3, upper=3.8> gamma;
}
model {
  // weakly informative priors
  mu ~ exponential(1/10);     // mean ~ 10 events
  gamma ~ uniform(1.3, 3.8);

  vector[N] mix = (mu / N) .* (Sspace .* Sen) + (1 - mu / N) .* Ben;
  target += sum(log(mix + 1e-300));
}
generated quantities {
  // Posterior predictive log-lik (for LOO if desired)
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = log( (mu/N) * Sspace[n] * Sen[n] + (1 - mu/N) * Ben[n] + 1e-300 );
}
"

run_bayes_single <- function(src, df, kde_bkg, iter=2000, chains=4, mode_energy="fallback_pl") {
  xmin <- min(df$Elog10); xmax <- max(df$Elog10)
  Ssp  <- compute_S_space(src$ra_deg, src$dec_deg, df$AngErr, df$RA, df$Dec, kernel="gauss")
  Sen  <- eval_S_energy(df$Elog10, df$Dec, gamma=2.2, kde_bkg, xmin, xmax, mode_energy) # initialized gamma is ignored in Stan; we pass vector but Stan reweights via gamma param? -> Trick: pass Sen at gamma=2.2 and re-evaluate in model? 
  # NOTE: Because Sen depends on gamma, we cannot pass a fixed Sen. We'll recompute Sen inside Stan via closed form.
  # For simplicity, we now pass Elog10 and bounds so Stan can compute Sen analytically.
  
  # Rebuild Stan with energy done inside:
  stan_code <- "
  data {
    int<lower=1> N;
    vector<lower=0>[N] Sspace;
    vector[N] elog10;
    vector[N] dec;
    real xmin;
    real xmax;
    vector<lower=0>[N] Ben;  // background energy density from KDE (fixed)
  }
  parameters {
    real<lower=0> mu;
    real<lower=1.3, upper=3.8> gamma;
  }
  transformed parameters {
    vector<lower=0>[N] Sen;
    {
      real a;
      a = 1 - gamma;
      if (fabs(a) < 1e-8) {
        // nearly flat in x
        for (n in 1:N) Sen[n] = 1/(xmax - xmin);
      } else {
        real Z = (pow(10, a*xmax) - pow(10, a*xmin)) / (a * log(10));
        for (n in 1:N) Sen[n] = pow(10, a*elog10[n]) / Z;
      }
    }
  }
  model {
    mu ~ exponential(0.1);
    gamma ~ uniform(1.3, 3.8);
    vector[N] mix = (mu / N) .* (Sspace .* Sen) + (1 - mu / N) .* Ben;
    target += sum(log(mix + 1e-300));
  }
  generated quantities {
    vector[N] log_lik;
    for (n in 1:N)
      log_lik[n] = log( (mu/N) * Sspace[n] * Sen[n] + (1 - mu/N) * Ben[n] + 1e-300 );
  }
  "
  
  # Background energy density fixed from KDE:
  Ben <- eval_B_energy(df$Elog10, df$Dec, kde_bkg)
  
  standata <- list(
    N = nrow(df),
    Sspace = as.vector(Ssp),
    elog10 = as.vector(df$Elog10),
    dec = as.vector(df$Dec),
    xmin = xmin, xmax = xmax,
    Ben = as.vector(Ben)
  )
  
  fit <- NULL
  if (!is.null(cmdstanr::cmdstan_version(error = FALSE))) {
    mod <- cmdstanr::cmdstan_model(write_stan_file(stan_code))
    fit <- mod$sample(data=standata, iter_warmup=1000, iter_sampling=iter, chains=chains,
                      refresh=200, adapt_delta=0.9, max_treedepth=12)
    draws <- fit$draws()
    summ  <- posterior::summarise_draws(draws, "mu","gamma")
    posterior::write_draws(fit$draws(), file.path(OUTDIR,"diagnostics", glue("draws_single_{str_replace_all(src$name,'[^A-Za-z0-9]+','_')}.nc")))
  } else {
    cli::cli_alert_warning("CmdStan not available. Falling back to rstan (slower).")
    rstan::rstan_options(auto_write = TRUE)
    sm <- rstan::stan_model(model_code = stan_code)
    fit <- rstan::sampling(sm, data = standata, iter = 2000, chains = chains, refresh = 100)
    summ <- rstan::summary(fit, pars = c("mu","gamma"))$summary %>% 
      as_tibble(rownames = "param")
  }
  
  # Plots
  if (!is.null(cmdstanr::cmdstan_version(error=FALSE))) {
    df_draws <- as_draws_df(fit$draws())
    p1 <- bayesplot::mcmc_areas(df_draws, pars = c("mu","gamma")) + ggtitle(glue("Posterior: {src$name}"))
  } else {
    df_draws <- as.data.frame(rstan::extract(fit, pars=c("mu","gamma")))
    p1 <- ggplot(df_draws, aes(mu)) + geom_density() + theme_paper +
      ggtitle(glue("Posterior mu: {src$name}"))
  }
  ggsave(file.path(OUTDIR,"figures", glue("posterior_{str_replace_all(src$name,'[^A-Za-z0-9]+','_')}.png")),
         p1, width=6.8, height=4, dpi=220)
  
  list(fit=fit, summary=summ)
}

# Run Bayesian fit for the top frequentist source by TS
best_src <- res_tbl %>% slice_max(TS, n=1) %>% slice(1)
cli::cli_h1("Bayesian fit (single source): {best_src$source}")
post <- run_bayes_single(best_src, df, bkg_kde, iter=2000, chains=4, mode_energy="fallback_pl")

# Save posterior summary
if (inherits(post$summary, "data.frame")) {
  write_csv(as_tibble(post$summary), file.path(OUTDIR,"tables","bayes_single_summary.csv"))
} else {
  posterior::write_table(post$summary, file.path(OUTDIR,"tables","bayes_single_summary.csv"))
}

# ---------- 13) Optional: Hierarchical stacking over catalog ----------
stan_code_hier <- "
data {
  int<lower=1> N;              // events
  int<lower=1> S;              // sources
  array[N] int<lower=1, upper=S> sidx;  // which source weight applies to event n
  vector<lower=0>[N] Sspace;   // already includes per-source angular separations
  vector<lower=0>[N] Ben;      
  vector[N] elog10;           
  real xmin;
  real xmax;
  vector<lower=0>[S] w;        // stacking weights (normalized or not)
}
parameters {
  real eta;                     // log-mean of mu_s
  real<lower=0> tau;            // log-sd of mu_s
  vector<lower=0>[S] mu_s;      // per-source expected counts
  real<lower=1.3, upper=3.8> gamma;
}
transformed parameters {
  vector<lower=0>[N] Sen;
  {
    real a = 1 - gamma;
    if (fabs(a) < 1e-8) {
      for (n in 1:N) Sen[n] = 1/(xmax - xmin);
    } else {
      real Z = (pow(10, a*xmax) - pow(10, a*xmin)) / (a * log(10));
      for (n in 1:N) Sen[n] = pow(10, a*elog10[n]) / Z;
    }
  }
}
model {
  tau ~ normal(0, 1);
  eta ~ normal(log(10), 1);   // prior centers around ~10 counts per source on average
  mu_s ~ lognormal(eta, tau);

  vector[N] mix;
  for (n in 1:N) {
    real mu_stack = dot_product(w, mu_s);              // total expected signal across catalog
    real frac     = mu_stack / N;
    mix[n] = frac * Sspace[n] * Sen[n] + (1 - frac) * Ben[n];
  }
  target += sum(log(mix + 1e-300));
}
"

run_bayes_hier <- function(cat_df, df, kde_bkg, w = NULL, iter=2000, chains=4) {
  S <- nrow(cat_df)
  if (is.null(w)) w <- rep(1/S, S)
  # Build Sspace stacked: sum_s w_s * Sspace_s evaluated per event, but we also pass sidx for clarity
  Sspace_stack <- rep(0, nrow(df))
  for (s in seq_len(S)) {
    Sspace_stack <- Sspace_stack + w[s] * compute_S_space(cat_df$ra_deg[s], cat_df$dec_deg[s],
                                                          df$AngErr, df$RA, df$Dec, kernel="gauss")
  }
  Ben <- eval_B_energy(df$Elog10, df$Dec, kde_bkg)
  standata <- list(
    N = nrow(df),
    S = S,
    sidx = rep(1, nrow(df)),   # not used in current parametrization
    Sspace = as.vector(Sspace_stack),
    Ben = as.vector(Ben),
    elog10 = as.vector(df$Elog10),
    xmin = min(df$Elog10), xmax = max(df$Elog10),
    w = as.vector(w)
  )
  if (!is.null(cmdstanr::cmdstan_version(error=FALSE))) {
    mod <- cmdstanr::cmdstan_model(write_stan_file(stan_code_hier))
    fit <- mod$sample(data=standata, iter_warmup=1000, iter_sampling=iter, chains=chains,
                      refresh=200, adapt_delta=0.9, max_treedepth=12)
    draws <- fit$draws()
    summ  <- posterior::summarise_draws(draws, dplyr::all_of(c("eta","tau","gamma")))
    posterior::write_draws(draws, file.path(OUTDIR,"diagnostics","draws_hier_catalog.nc"))
  } else {
    cli::cli_alert_warning("CmdStan not available; skipping hierarchical Bayesian stacking.")
    return(NULL)
  }
  list(fit=fit, summary=summ)
}

cli::cli_h1("Hierarchical stacking (equal weights)")
hier <- run_bayes_hier(cat_df, df, bkg_kde, w = rep(1/nrow(cat_df), nrow(cat_df)), iter=1500, chains=4)
if (!is.null(hier)) {
  posterior::write_table(hier$summary, file.path(OUTDIR,"tables","bayes_hier_summary.csv"))
}

cli::cli_h1("Done.")
cli::cli_alert_success("Figures: {file.path(OUTDIR,'figures')}")
cli::cli_alert_success("Tables:  {file.path(OUTDIR,'tables')}")
cli::cli_alert_success("Diag:     {file.path(OUTDIR,'diagnostics')}")

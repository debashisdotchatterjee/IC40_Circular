# icecube_allsky_hotspots.R
# ALL-SKY hotspot search using ONLY the IC40 event dataset (no external catalog).
# Robust column detection; beautiful plots; TS map; RA-scrambled calibration; Bayesian HMC.

# -------------------- 0) Packages & Setup --------------------
req <- c(
  "data.table","dplyr","ggplot2","ggpointdensity","scales","readr","stringr",
  "lubridate","purrr","tidyr","glue","fs","cli","progress","pbapply",
  "matrixStats","viridis","circular","knitr","kableExtra","posterior",
  "cmdstanr","rstan","hexbin"
)
miss <- setdiff(req, rownames(installed.packages()))
if(length(miss)) install.packages(miss, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))
options(mc.cores = parallel::detectCores())

# Output folders
stamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- glue("allsky_outputs_{stamp}")
for (d in c("figures","tables","diagnostics","cache")) fs::dir_create(file.path(OUTDIR,d))
cli::cli_h1("Outputs -> {OUTDIR}")

# -------------------- 1) Helpers --------------------
deg2rad <- function(x) x*pi/180; rad2deg <- function(x) x*180/pi
wrap360 <- function(x) (x %% 360)
angsep_deg <- function(ra1, dec1, ra2, dec2){
  ra1r <- deg2rad(ra1); dec1r <- deg2rad(dec1)
  ra2r <- deg2rad(ra2); dec2r <- deg2rad(dec2)
  d <- acos(pmin(1, pmax(-1, sin(dec1r)*sin(dec2r)+cos(dec1r)*cos(dec2r)*cos(ra1r-ra2r))))
  rad2deg(d)
}
S_space_gauss <- function(dpsi_deg, sigma_deg){
  s <- pmax(sigma_deg, 0.2)
  (1/(2*pi*s^2)) * exp(-(dpsi_deg^2)/(2*s^2))
}
S_energy_pl_log10 <- function(x, gamma, xmin, xmax){
  a <- 1-gamma
  if (abs(a)<1e-6) return(rep(1/(xmax-xmin), length(x)))
  Z <- (10^(a*xmax) - 10^(a*xmin)) / (a*log(10))
  10^(a*x)/Z
}
gmst_from_mjd <- function(mjd){ JD <- mjd + 2400000.5; T <- (JD-2451545.0)/36525
gmst_sec <- 67310.54841 + (876600*3600+8640184.812866)*T + 0.093104*T^2 - 6.2e-6*T^3
wrap360(gmst_sec/240)
}
ra_dec_from_azzen_pole <- function(az_deg, zen_deg, mjd){
  dec <- -(90 - zen_deg); H <- wrap360(az_deg); lst <- gmst_from_mjd(mjd)
  ra  <- wrap360(lst - H); tibble(RA=ra, Dec=dec)
}

# -------------------- 2) Load events (robust column detection) --------------------
# Try local, then /mnt/data
infile <- if (file.exists("IC40_exp.csv")) "IC40_exp.csv" else "/mnt/data/IC40_exp.csv"
stopifnot(file.exists(infile))
raw <- suppressWarnings(data.table::fread(infile))
names(raw) <- stringr::str_trim(names(raw))

# Standardize headers
std <- function(nm){
  nm <- str_replace_all(nm, "\\[.*?\\]","")    # strip bracketed units
  nm <- str_replace_all(nm, "[^A-Za-z0-9]+","_")
  tolower(nm)
}
names(raw) <- std(names(raw))

# Helper: find first name matching ANY of the regexes (in order)
find_by_patterns <- function(nms, patterns){
  for (p in patterns) {
    hits <- which(grepl(p, nms, perl=TRUE))
    if (length(hits)) return(nms[hits[1]])
  }
  NA_character_
}

nms <- names(raw)

# Detect columns
col_mjd <- find_by_patterns(nms, c("^mjd(_days)?$","mjd"))
col_elg <- find_by_patterns(nms, c("elog10","log10.*e.*gev.*","^log10e.*","^log_?10.*e.*gev"))
col_ang <- find_by_patterns(nms, c("^angerr(_deg)?$","ang.*err"))
col_ra  <- find_by_patterns(nms, c("^ra(_deg)?$","^ra$","right_?asc"))
col_dc  <- find_by_patterns(nms, c("^dec(_deg)?$","^dec$","decl"))
col_az  <- find_by_patterns(nms, c("^azimuth(_deg)?$","^azimuth$","^az$"))
col_zn  <- find_by_patterns(nms, c("^zenith(_deg)?$","^zenith$","^zen$"))

needed <- c(col_mjd,col_elg,col_ang,col_ra,col_dc)
if (any(is.na(needed))) {
  cli::cli_abort(c(
    "x" = "Could not find required columns after header cleaning.",
    "i" = glue("Detected -> MJD:{col_mjd}  Elog10:{col_elg}  AngErr:{col_ang}  RA:{col_ra}  Dec:{col_dc}"),
    "i" = "Please open the CSV and check the column names."
  ))
}

# Build main DF (Az/Zen optional)
df <- raw %>% transmute(
  MJD   = as.numeric(.data[[col_mjd]]),
  Elog10= as.numeric(.data[[col_elg]]),
  AngErr= as.numeric(.data[[col_ang]]),
  RA    = wrap360(as.numeric(.data[[col_ra]])),
  Dec   = as.numeric(.data[[col_dc]]),
  Az    = if (!is.na(col_az)) as.numeric(.data[[col_az]]) else NA_real_,
  Zen   = if (!is.na(col_zn)) as.numeric(.data[[col_zn]]) else NA_real_
) %>% filter(is.finite(MJD), is.finite(Elog10), is.finite(AngErr), is.finite(RA), is.finite(Dec))

cli::cli_alert_success("Loaded {nrow(df)} events.")
cli::cli_inform(glue("Columns -> MJD:{col_mjd}  Elog10:{col_elg}  AngErr:{col_ang}  RA:{col_ra}  Dec:{col_dc}  Az:{col_az %||% 'NA'}  Zen:{col_zn %||% 'NA'}"))

# -------------------- 3) Descriptives & pretty plots --------------------
desc <- tibble(
  MJD_min=min(df$MJD), MJD_max=max(df$MJD),
  E_med=median(df$Elog10), E_q05=quantile(df$Elog10,.05), E_q95=quantile(df$Elog10,.95),
  S_med=median(df$AngErr), S_q05=quantile(df$AngErr,.05), S_q95=quantile(df$AngErr,.95),
  frac_sigma_lt1=mean(df$AngErr<1)
)
readr::write_csv(desc, file.path(OUTDIR,"tables","descriptives.csv"))

theme_paper <- theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
        legend.position="right", strip.text=element_text(face="bold"))

# E distribution
pE <- ggplot(df,aes(Elog10))+
  geom_histogram(aes(y=after_stat(density)),bins=80,alpha=.85)+
  geom_density(linewidth=1)+
  labs(title="Reconstructed log10(E/GeV)",x=expression(log[10](E/GeV)),y="Density") + theme_paper
ggsave(file.path(OUTDIR,"figures","elog10_hist.png"),pE,width=7.4,height=4.4,dpi=220)

# AngErr distribution
pS <- ggplot(df,aes(AngErr))+
  geom_histogram(bins=80,alpha=.85)+
  labs(title="Per-event angular uncertainty",x="AngErr [deg]",y="Count")+theme_paper
ggsave(file.path(OUTDIR,"figures","angerr_hist.png"),pS,width=7.4,height=4.4,dpi=220)

# Sky density
pSky <- ggplot(df, aes(RA, Dec))+
  ggpointdensity::geom_pointdensity(size=0.25)+
  scale_color_viridis_c(option="plasma")+
  coord_cartesian(xlim=c(0,360), ylim=c(-90,90), expand=FALSE)+
  labs(title="IC40 sky map (event density)", x="RA [deg]", y="Dec [deg]", color="Local\nDensity")+
  theme_paper
ggsave(file.path(OUTDIR,"figures","sky_density.png"), pSky, width=8.2, height=4.6, dpi=220)

# Energy vs AngErr hexbin
pHex <- ggplot(df, aes(Elog10, AngErr))+
  geom_hex(bins=60)+
  scale_fill_viridis(option="turbo")+
  labs(title="Energy vs Angular Error", x=expression(log[10](E/GeV)), y="AngErr [deg]", fill="Count")+
  theme_paper
ggsave(file.path(OUTDIR,"figures","energy_vs_angerr_hex.png"), pHex, width=7.4, height=4.6, dpi=220)

# RA uniformity
ra_circ <- circular::circular(deg2rad(df$RA), type="angles", units="radians", template="geographics")
ra_uni  <- tibble(test="Rayleigh_RA_uniformity", p_value=as.numeric(circular::rayleigh.test(ra_circ)$p.value))
write_csv(ra_uni, file.path(OUTDIR,"tables","ra_uniformity.csv"))

# -------------------- 4) Background energy KDE (declination bands) --------------------
make_bkg_kde <- function(df, nbands=10){
  qs <- seq(0,1,length.out=nbands+1); cuts <- as.numeric(quantile(df$Dec,qs,na.rm=TRUE))
  bands <- cut(df$Dec, breaks=unique(cuts), include.lowest=TRUE)
  mods <- vector("list", length(levels(bands)))
  for (k in seq_along(mods)){
    x <- df$Elog10[bands==levels(bands)[k]]
    mods[[k]] <- density(x, n=2048, bw="nrd0", from=min(x), to=max(x))
  }
  list(cuts=unique(cuts), mods=mods)
}
eval_bkg_kde <- function(elog10, dec, kde){
  b <- pmin(findInterval(dec, kde$cuts, all.inside=TRUE), length(kde$mods))
  out <- numeric(length(elog10))
  for (k in seq_along(kde$mods)){
    idx <- which(b==k)
    if (length(idx)) out[idx] <- approx(kde$mods[[k]]$x, kde$mods[[k]]$y, elog10[idx], rule=2)$y
  }
  pmax(out, .Machine$double.eps)
}
cli::cli_h1("Fitting background energy KDE …")
bkg_kde <- make_bkg_kde(df, nbands=10); saveRDS(bkg_kde, file.path(OUTDIR,"cache","bkg_kde.rds"))

# -------------------- 5) Energy-weighted pre-scan --------------------
gamma0 <- 2.2; xmin <- min(df$Elog10); xmax <- max(df$Elog10)
S_en0  <- S_energy_pl_log10(df$Elog10, gamma0, xmin, xmax)
B_en   <- eval_bkg_kde(df$Elog10, df$Dec, bkg_kde)
w_en   <- pmax(S_en0 / B_en, 1e-6)

# Grid resolution (tune)
GRID_RA_STEP  <- 1.0    # deg
GRID_DEC_STEP <- 1.0    # deg
trial_RA  <- seq(0, 360 - GRID_RA_STEP, by=GRID_RA_STEP)
trial_Dec <- seq(-90, 90, by=GRID_DEC_STEP)
grid <- expand.grid(ra_deg=trial_RA, dec_deg=trial_Dec)
cli::cli_alert_info("Pre-scan grid -> {nrow(grid)} directions.")

sigma0 <- as.numeric(quantile(df$AngErr, 0.5))  # median AngErr
RADIUS <- 3.0 * sigma0

cli::cli_h1("Pre-scan energy-weighted map")
pb <- progress::progress_bar$new(total=nrow(grid),
                                 format="  pre-scan [:bar] :percent :current/:total eta: :eta", clear=FALSE, width=70)

pre_scan_vals <- pbsapply(seq_len(nrow(grid)), function(i){
  ra0  <- grid$ra_deg[i]; dec0 <- grid$dec_deg[i]
  dec_sel <- abs(df$Dec - dec0) <= (RADIUS + 1.0)
  dra <- pmin(abs(df$RA - ra0), 360 - abs(df$RA - ra0))
  ra_sel <- dra <= (RADIUS + 1.0)
  idx <- which(dec_sel & ra_sel)
  if (!length(idx)) { pb$tick(); return(0) }
  dpsi <- angsep_deg(df$RA[idx], df$Dec[idx], ra0, dec0)
  idx2 <- idx[ dpsi <= (RADIUS + 2.0) ]
  if (!length(idx2)) { pb$tick(); return(0) }
  dens <- S_space_gauss(angsep_deg(df$RA[idx2], df$Dec[idx2], ra0, dec0), sigma0)
  val  <- sum(w_en[idx2] * dens)
  pb$tick(); val
})
grid$pre_score <- as.numeric(pre_scan_vals)

pPre <- ggplot(grid, aes(ra_deg, dec_deg, fill=pre_score)) +
  geom_raster(interpolate=TRUE) + scale_fill_viridis(option="plasma") +
  coord_cartesian(xlim=c(0,360), ylim=c(-90,90), expand=FALSE) +
  labs(title="All-sky energy-weighted pre-scan", x="RA [deg]", y="Dec [deg]", fill="score") +
  theme_paper
ggsave(file.path(OUTDIR,"figures","allsky_prescan_heatmap.png"), pPre, width=9, height=4.8, dpi=220)

# -------------------- 6) Select top-K hotspots --------------------
TOP_K <- 50
top_idx <- order(grid$pre_score, decreasing=TRUE)[1:TOP_K]
cands  <- grid[top_idx,] %>% mutate(name = glue("hotspot_{row_number()}"))
readr::write_csv(cands, file.path(OUTDIR,"tables","hotspot_candidates.csv"))

# -------------------- 7) Full unbinned MLE at each hotspot --------------------
eval_B_space <- function(dec) rep(1, length(dec)) # uniform in RA (time-integrated)
negloglik_pt <- function(par, ra0, dec0, df, kde){
  mu    <- pmax(par[1], 0)
  gamma <- max(1.3, min(par[2], 3.8))
  N     <- nrow(df)
  dpsi  <- angsep_deg(df$RA, df$Dec, ra0, dec0)
  Ssp   <- S_space_gauss(dpsi, df$AngErr)
  Sen   <- S_energy_pl_log10(df$Elog10, gamma, xmin=min(df$Elog10), xmax=max(df$Elog10))
  Ben   <- eval_bkg_kde(df$Elog10, df$Dec, kde)
  Bsp   <- eval_B_space(df$Dec)
  S_i   <- Ssp * Sen
  B_i   <- Bsp * Ben
  Li    <- (mu/N) * S_i + (1 - mu/N) * B_i
  -sum(log(pmax(Li, 1e-300)))
}
fit_hotspot <- function(ra0, dec0, df, kde){
  par0 <- c(mu=10, gamma=2.2)
  fit  <- nlminb(par0, negloglik_pt, lower=c(0,1.2), upper=c(Inf,4.2),
                 ra0=ra0, dec0=dec0, df=df, kde=kde)
  nll_hat <- fit$objective
  nll0    <- negloglik_pt(c(0,2.2), ra0=ra0, dec0=dec0, df=df, kde=kde)
  TS <- 2*(nll0 - nll_hat)
  list(mu_hat=max(fit$par[1],0),
       gamma_hat=max(1.3,min(fit$par[2],3.8)),
       TS=max(TS,0),
       conv=(fit$convergence==0))
}

cli::cli_h1("Refining top-{TOP_K} hotspots (MLE)")
MLE_res <- vector("list", nrow(cands))
pb2 <- progress::progress_bar$new(total=nrow(cands),
                                  format="  refine [:bar] :percent :current/:total eta: :eta", width=70)
for (i in seq_len(nrow(cands))){
  hs <- cands[i,]
  fr <- fit_hotspot(hs$ra_deg, hs$dec_deg, df, bkg_kde)
  MLE_res[[i]] <- tibble(
    name=hs$name, ra=hs$ra_deg, dec=hs$dec_deg,
    mu_hat=fr$mu_hat, gamma_hat=fr$gamma_hat, TS=fr$TS, converged=fr$conv
  )
  pb2$tick()
}
MLE_tbl <- bind_rows(MLE_res) %>% arrange(desc(TS))
write_csv(MLE_tbl, file.path(OUTDIR,"tables","hotspots_mle_results.csv"))

# -------------------- 8) Post-trial calibration (RA-scrambling) --------------------
scramble_ra_uniform <- function(n) runif(n, 0, 360)
scramble_time_recompute_radec <- function(df){
  if (all(!is.finite(df$Az)) || all(!is.finite(df$Zen))) return(scramble_ra_uniform(nrow(df)))
  mjd_star <- sample(df$MJD, size=nrow(df), replace=TRUE)
  ra_star  <- ra_dec_from_azzen_pole(df$Az, df$Zen, mjd_star)$RA
  if (any(!is.finite(ra_star))) ra_star <- scramble_ra_uniform(nrow(df))
  ra_star
}
calibrate_one <- function(ra0, dec0, df, kde, n_sims=500, method=c("uniformRA","time")){
  method <- match.arg(method)
  N <- nrow(df); Bsp <- eval_B_space(df$Dec)
  xmin <- min(df$Elog10); xmax <- max(df$Elog10)
  Ben <- eval_bkg_kde(df$Elog10, df$Dec, kde)
  pb <- progress::progress_bar$new(total=n_sims,
                                   format=glue("  scramble @RA={round(ra0,2)},Dec={round(dec0,2)} [:bar] :percent eta: :eta"),
                                   clear=FALSE, width=70)
  TSnull <- pbsapply(seq_len(n_sims), function(b){
    RA_star <- if (method=="uniformRA") scramble_ra_uniform(N) else scramble_time_recompute_radec(df)
    dpsi <- angsep_deg(RA_star, df$Dec, ra0, dec0)
    Ssp  <- S_space_gauss(dpsi, df$AngErr)
    f <- nlminb(c(0,2.2),
                function(par){
                  mu    <- pmax(par[1], 0); gamma <- max(1.3,min(par[2],3.8))
                  Sen   <- S_energy_pl_log10(df$Elog10, gamma, xmin, xmax)
                  S_i   <- Ssp*Sen; B_i <- Bsp*Ben
                  Li    <- (mu/N)*S_i + (1 - mu/N)*B_i
                  -sum(log(pmax(Li,1e-300)))
                },
                lower=c(0,1.2), upper=c(Inf,4.2))
    nll0 <- {
      Sen0 <- S_energy_pl_log10(df$Elog10, 2.2, xmin, xmax)
      S_i0 <- Ssp*Sen0; B_i0 <- Bsp*Ben
      -sum(log(pmax(B_i0, 1e-300)))
    }
    TS <- 2*(nll0 - f$objective)
    pb$tick(); max(TS,0)
  })
  TSnull
}

N_CALIB <- 500
TOP_CAL <- min(10, nrow(MLE_tbl))

cli::cli_h1("Post-trial calibration for top {TOP_CAL}")
cal_list <- vector("list", TOP_CAL)
for (i in seq_len(TOP_CAL)){
  hs <- MLE_tbl[i,]
  cal_list[[i]] <- calibrate_one(hs$ra, hs$dec, df, bkg_kde, n_sims=N_CALIB,
                                 method=if (any(is.finite(df$Az)&is.finite(df$Zen))) "time" else "uniformRA")
  # Save null hist
  dnull <- tibble(TS = cal_list[[i]])
  pnull <- ggplot(dnull, aes(TS))+geom_histogram(bins=60, fill="grey70", color="white")+
    geom_vline(xintercept=hs$TS, color="red", linewidth=1)+
    labs(title=glue("{hs$name}: RA-scrambled null TS"), x="TS", y="Count")+theme_paper
  ggsave(file.path(OUTDIR,"figures", glue("TSnull_{hs$name}.png")), pnull, width=7.2, height=4.4, dpi=220)
}
MLE_tbl$p_post <- NA_real_
MLE_tbl$p_post[seq_len(TOP_CAL)] <- sapply(seq_len(TOP_CAL), function(i) mean(cal_list[[i]] >= MLE_tbl$TS[i]))
write_csv(MLE_tbl, file.path(OUTDIR,"tables","hotspots_mle_with_p.csv"))

# -------------------- 9) TS map & summary plots --------------------
# Put TS on a narrow grid around candidates for a pretty map:
grid$TS <- NA_real_
# Match candidate rows to grid rows (exact positions)
key <- paste(cands$ra_deg, cands$dec_deg)
idx <- match(paste(grid$ra_deg, grid$dec_deg), key)
grid$TS[!is.na(idx)] <- MLE_tbl$TS[idx[!is.na(idx)]]

pTS <- ggplot(grid, aes(ra_deg, dec_deg, fill=TS)) +
  geom_raster(na.rm=TRUE, interpolate=TRUE) +
  scale_fill_viridis(option="magma", na.value="grey90") +
  coord_cartesian(xlim=c(0,360), ylim=c(-90,90), expand=FALSE) +
  labs(title="Refined TS at top pre-scan peaks", x="RA [deg]", y="Dec [deg]", fill="TS") +
  theme_paper
ggsave(file.path(OUTDIR,"figures","ts_refined_heatmap.png"), pTS, width=9, height=4.8, dpi=220)

best_tbl <- MLE_tbl %>%
  mutate(across(c(mu_hat,gamma_hat,TS,p_post), ~round(.,6))) %>%
  slice_head(n=TOP_CAL)
kt <- best_tbl %>% knitr::kable("html", caption="Top hotspots with post-trial p for top 10") %>%
  kableExtra::kable_styling(full_width=FALSE, bootstrap_options=c("striped","hover","condensed"))
htmltools::save_html(kt, file.path(OUTDIR,"tables","top_hotspots.html"))

# -------------------- 10) Bayesian HMC on top few --------------------
stan_code_single <- "
data {
  int<lower=1> N;
  vector<lower=0>[N] Sspace;
  vector[N] elog10;
  real xmin;
  real xmax;
  vector<lower=0>[N] Ben;
}
parameters {
  real<lower=0> mu;
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
  mu ~ exponential(0.1);
  gamma ~ uniform(1.3, 3.8);
  vector[N] mix = (mu/N) .* (Sspace .* Sen) + (1 - mu/N) .* Ben;
  target += sum(log(mix + 1e-300));
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = log( (mu/N) * Sspace[n] * Sen[n] + (1 - mu/N) * Ben[n] + 1e-300 );
}
"

run_bayes_pos <- function(ra0, dec0, df, kde, iter=2000, chains=4){
  Ben <- eval_bkg_kde(df$Elog10, df$Dec, kde)
  dpsi <- angsep_deg(df$RA, df$Dec, ra0, dec0)
  Ssp  <- S_space_gauss(dpsi, df$AngErr)
  standata <- list(
    N = nrow(df), Sspace = as.vector(Ssp),
    elog10 = as.vector(df$Elog10),
    xmin = min(df$Elog10), xmax = max(df$Elog10),
    Ben = as.vector(Ben)
  )
  if (!is.null(cmdstanr::cmdstan_version(error=FALSE))) {
    mod <- cmdstanr::cmdstan_model(write_stan_file(stan_code_single))
    fit <- mod$sample(data=standata, iter_warmup=1000, iter_sampling=iter, chains=chains,
                      refresh=200, adapt_delta=0.9, max_treedepth=12)
    draws <- fit$draws()
    smry  <- posterior::summarise_draws(draws, c("mu","gamma"))
    return(list(fit=fit, summary=smry, draws=draws))
  } else {
    rstan::rstan_options(auto_write=TRUE)
    sm <- rstan::stan_model(model_code=stan_code_single)
    fit <- rstan::sampling(sm, data=standata, iter=iter, chains=chains, refresh=100)
    smry <- as_tibble(rstan::summary(fit, pars=c("mu","gamma"))$summary, rownames="param")
    return(list(fit=fit, summary=smry, draws=NULL))
  }
}

TOP_BAYES <- min(5, nrow(MLE_tbl))
cli::cli_h1("Bayesian HMC for top {TOP_BAYES} hotspots")
bayes_list <- vector("list", TOP_BAYES)
for (i in seq_len(TOP_BAYES)){
  hs <- MLE_tbl[i,]
  bay <- run_bayes_pos(hs$ra, hs$dec, df, bkg_kde, iter=2000, chains=4)
  bayes_list[[i]] <- list(name=hs$name, ra=hs$ra, dec=hs$dec, summary=bay$summary)
  # posterior plots if cmdstanr present
  if (!is.null(cmdstanr::cmdstan_version(error=FALSE)) && !is.null(bay$draws)) {
    dr <- as_draws_df(bay$draws)
    p_mu <- bayesplot::mcmc_areas(dr, pars="mu") + ggtitle(glue("{hs$name}: posterior mu"))
    p_ga <- bayesplot::mcmc_areas(dr, pars="gamma") + ggtitle(glue("{hs$name}: posterior gamma"))
    ggsave(file.path(OUTDIR,"figures", glue("post_mu_{hs$name}.png")), p_mu, width=6.5, height=3.6, dpi=220)
    ggsave(file.path(OUTDIR,"figures", glue("post_gamma_{hs$name}.png")), p_ga, width=6.5, height=3.6, dpi=220)
  }
}
bx <- purrr::imap_dfr(bayes_list, ~{ sm <- as_tibble(.x$summary); sm$name <- .x$name; sm$ra <- .x$ra; sm$dec <- .x$dec; sm })
write_csv(bx, file.path(OUTDIR,"tables","bayes_top_summaries.csv"))

cli::cli_alert_success("Done. Figures: {file.path(OUTDIR,'figures')}  Tables: {file.path(OUTDIR,'tables')}")

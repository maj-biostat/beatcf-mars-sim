suppressPackageStartupMessages(library("cmdstanr"))
suppressPackageStartupMessages(library("posterior"))
suppressPackageStartupMessages(library("config"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("ggthemes"))
suppressPackageStartupMessages(library("qs"))
suppressPackageStartupMessages(library("extraDistr"))
suppressPackageStartupMessages(library("pbapply"))
suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library(mvtnorm))

ix <- 1

stan_model_code <- "
data {
  int<lower=1> N;                 // observations
  int<lower=1> P; // matrix cols
  
  matrix[N, P] X;
  array[N] int y; // 1 none, 2 mild, 3 severe

  int ix_trt_2;
  int  ix_trt_3 ;
  int  ix_prev_2;
  int  ix_prev_3;
  int  ix_time_1;
  int  ix_time_2;
  int  ix_prev_time_2;
  int  ix_prev_time_3;
  int  ix_trt_time_2;
  int  ix_trt_time_3;
  
  real mu_days;
  real sd_days;
  
}
transformed data{
  vector[1] zero1 = rep_vector(0, 1);
}
parameters {
  ordered[2] alpha;
  vector[P] b;
}
transformed parameters {
}
model {

  // priors
  alpha ~ normal(0, 2);
  b ~ normal(0, 1.5);

  // likelihood
  
  vector[N] eta = X * b;
  y ~ ordered_logistic(eta, alpha);

}
generated quantities{
  
  // see sim18-back-transform-note.qmd
    vector[2] a;
  for (k in 1:2) {
    a[k] = alpha[k] - ((b[ix_time_2] * pow(mu_days,2)) / pow(sd_days, 2) - (b[ix_time_1] * mu_days) / sd_days);
  }

  vector[3] b_trt;
  b_trt[1] = 0.0;
  b_trt[2] = b[ix_trt_2] - (b[ix_trt_time_2] * mu_days) / sd_days;
  b_trt[3] = b[ix_trt_3] - (b[ix_trt_time_3] * mu_days) / sd_days;
  
  vector[3] b_prev;
  b_prev[1] = 0.0;
  b_prev[2] = b[ix_prev_2] - (b[ix_prev_time_2] * mu_days) / sd_days;
  b_prev[3] = b[ix_prev_3] - (b[ix_prev_time_3] * mu_days) / sd_days;
  
  real b_time_1 = (b[ix_time_1] / sd_days) - (2.0 * mu_days * b[ix_time_2]) / pow(sd_days, 2);
  real b_time_2 = b[ix_time_2] / pow(sd_days, 2);
  
  vector[3] b_prev_time;
  b_prev_time[1] = 0.0;
  b_prev_time[2] = b[ix_prev_time_2] / sd_days;
  b_prev_time[3] = b[ix_prev_time_3] / sd_days;
  
  vector[3] b_trt_time;
  b_trt_time[1] = 0.0;
  b_trt_time[2] = b[ix_trt_time_2] / sd_days;
  b_trt_time[3] = b[ix_trt_time_3] / sd_days;
  
}
"

stan_file <- write_stan_file(stan_model_code)
m_1 <- cmdstanr::cmdstan_model(stan_file)

output_dir_mcmc <- paste0(getwd(), "/tmp")

sim18_cohort <- function(l_spec){
  
  id_cohort <- l_spec$is:l_spec$ie
  N_cohort  <- length(id_cohort)
  n_day     <- l_spec$max_day + 1
  
  ## subject-level quantities
  trt_subj <- sample(
    l_spec$trt_lab[l_spec$trt_active],
    N_cohort,
    replace = TRUE
  )
  
  ## cohort data
  d_cohort <- data.table(
    id  = rep(id_cohort, each = n_day),
    day = rep(0:l_spec$max_day, times = N_cohort),
    trt = rep(trt_subj, each = n_day)
  )
  
  d_cohort[, t_0 := l_spec$t_0[id]]
  
  ## work directly on vectors
  state <- integer(nrow(d_cohort))
  
  ## initial states
  state[seq(1L, by = n_day, length.out = N_cohort)] <-
    sample(
      l_spec$state_opts,
      N_cohort,
      replace = TRUE,
      prob = l_spec$p_init
    )
  
  ## cache coefficients
  b_trt        <- l_spec$b_trt
  b_prev       <- l_spec$b_prev
  b_prev_time  <- l_spec$b_prev_time
  b_trt_time <- l_spec$b_trt_time
  alpha       <- l_spec$alpha
  
  ## precompute time effect
  day_vals <- 0:l_spec$max_day
  time_eff <- l_spec$b_time_1 * day_vals + l_spec$b_time_2 * day_vals^2
  i <- 1
  
  for(i in seq_len(N_cohort)) {
    
    first <- (i - 1L) * n_day + 1L
    last  <- first + n_day - 1L
    
    trt_i <- trt_subj[i]
    
    for(r in (first + 1L):last) {
      
      prev <- state[r - 1L]
      tt   <- d_cohort$day[r]
      
      lp <-
        b_trt[trt_i] +
        b_prev[prev] +
        time_eff[tt + 1L] +
        b_prev_time[prev] * tt + 
        b_trt_time[trt_i] * tt 
      
      p0 <- plogis(alpha[1] - lp)
      p1 <- plogis(alpha[2] - lp)
      
      u <- runif(1)
      
      state[r] <-
        if (u < p0) {
          1L
        } else if (u < p1) {
          2L
        } else {
          3L
        }
    }
  }
  
  d_cohort[, state := state]
  
  d_obs <- copy(d_cohort)
  setorder(d_obs, id, day)
  
  # Need to compute previous state and previous day before we subset
  # otherwise we will end up referencing day 7 as the previous day from day 14
  d_obs[, `:=`(
    prev_state = data.table::shift(state, 1L),
    prev_day   = data.table::shift(day, 1L)
  ), by = id]
  
  d_obs <- d_obs[day %in% l_spec$visit_days]
  
  list(
    d_cohort = d_cohort,
    d_obs    = d_obs
  )
  
}


sim18_stan_data <- function(dd, l_spec){
  
  # For the day zero of onset, we have no prev state. We can either assume that they
  # were well the day before or just drop that observation and include it in day 1
  # of follow up. I do the latter.
  dd[, gap_len := day - prev_day]
  dd <- dd[day != 0]
  # gives 16 days of follow up per pt 14 whole days and then 1 obs in wk 3 and 4
  dd[, trt_idx := match(trt, l_spec$trt_lab)]
  
  dd[, x_time := copy(day)]
  # scaling makes this sample a lot faster but is a pain in the arse
  # for linear predictors with interactions
  # additionally means you need to back scale intercepts
  dd[, x_time := scale(x_time)]
  
  dd[, x_trt := factor(trt_idx)]
  dd[, x_prev := factor(prev_state, levels = l_spec$state_opts)]
  
  X <- model.matrix(~ x_trt + x_prev +
                      x_time + I(x_time^2) +
                      x_prev * x_time +
                      x_trt * x_time, 
                    data = dd)
  X_mod <- X[, -1]
  
  ld <- list(
    N  = nrow(dd),
    P = ncol(X_mod),
    X = X_mod,
    y = dd$state,
    # Indexing parameters within design matrix.
    # First list any variable with no dependency on time
    ix_trt_2 = 1,
    ix_trt_3 = 2,
    ix_prev_2 = 3,
    ix_prev_3 = 4,
    ix_time_1 = 5,
    ix_time_2 = 6,
    ix_prev_time_2 = 7,
    ix_prev_time_3 = 8,
    ix_trt_time_2 = 9,
    ix_trt_time_3 = 10,
    mu_days = mean(dd$day),
    sd_days = sd(dd$day)
  )
  stopifnot(names(X_mod)[ld$ix_trt_2] == "x_trt2")
  stopifnot(names(X_mod)[ld$ix_trt_3] == "x_trt3")
  stopifnot(names(X_mod)[ld$ix_prev_2] == "x_prev2")
  stopifnot(names(X_mod)[ld$ix_prev_3] == "x_prev3")
  stopifnot(names(X_mod)[ld$ix_time_1] == "x_time")
  stopifnot(names(X_mod)[ld$ix_time_2] == "I(x_time^2)")
  stopifnot(names(X_mod)[ld$ix_prev_time_2] == "x_prev2:x_time")
  stopifnot(names(X_mod)[ld$ix_prev_time_3] == "x_prev3:x_time")
  stopifnot(names(X_mod)[ld$ix_trt_time_2] == "x_trt2:x_time")
  stopifnot(names(X_mod)[ld$ix_trt_time_3] == "x_trt3:x_time")
  
  ld
}



sim18_transition_matrix <- function(day, trt, l_spec)
{
  
  P <- matrix(0,3,3)
  
  for(prev in l_spec$state_opts){
    
    lp <-
      l_spec$b_trt[trt] +
      l_spec$b_prev[prev] +
      l_spec$b_time_1 * day +
      l_spec$b_time_2 * day^2 +
      l_spec$b_prev_time[prev] * day +
      l_spec$b_trt_time[trt] * day 
    
    p1 <- plogis(l_spec$alpha[1] - lp)
    p2 <- plogis(l_spec$alpha[2] - lp) - p1
    p3 <- 1 - p1 - p2
    
    P[prev,] <- c(p1,p2,p3)
    
  }
  
  # from 
  rownames(P) <- l_spec$state_lab
  # to 
  colnames(P) <- l_spec$state_lab
  P
  
}


sim18_sop <- function(days = 1:28, l_spec)
{
  
  # has to start at day 1
  stopifnot(days[1] == 1)
  
  out <- matrix(NA, length(days) + 1, 3)
  out[1,] <- l_spec$p_init
  d_sop <- data.table()
  
  for(trt in l_spec$trt_lab){
    # starting point
    pi <- l_spec$p_init
    
    for(i in seq_along(days)){
      
      day = days[i]
      Pt <- sim18_transition_matrix(day, trt, l_spec)
      pi <- drop(pi %*% Pt)
      
      out[i + 1,] <- pi
    }
    
    d_sop <- rbind(
      d_sop, 
      data.table(
        trt = trt,
        day = c(0, days),
        none = out[,1],
        mild = out[,2],
        severe = out[,3]
      )
    )
  }
  
  d_sop
  
}


update_sim18_cfg <- function(l_spec){
  
  l_spec$trt_lab <- unlist(l_spec$trt_lab)
  l_spec$state_lab <- unlist(l_spec$state_lab)
  
  l_spec$alpha <- unlist(l_spec$alpha)
  l_spec$b_trt <- unlist(l_spec$b_trt)
  names(l_spec$b_trt) <- l_spec$trt_lab
  
  # *** Has to be converted to logical otherwise you are just going to be 
  # indexing trt 1
  l_spec$trt_active <- as.logical(l_spec$trt_active)
  names(l_spec$trt_active) <- l_spec$trt_lab
  
  l_spec$state_opts <- seq_along(l_spec$state_lab)
  names(l_spec$state_opts) <- l_spec$state_lab
  
  l_spec$visit_days <- unlist(l_spec$visit_days)
  l_spec$max_day <- max(l_spec$visit_days)
  
  l_spec$b_prev <- unlist(l_spec$b_prev)
  
  l_spec$b_trt_time <- unlist(l_spec$b_trt_time)
  names(l_spec$b_trt_time) <- l_spec$trt_lab
  
  l_spec$b_prev_time <- unlist(l_spec$b_prev_time)
  
  l_spec$p_init <- unlist(l_spec$p_init)
  
  l_spec$smry_pars <- c(
    "a", "b_trt", "b_prev", "b_time_1", "b_time_2", 
    "b_prev_time", 
    "b_trt_time")
  
  l_spec$full_pars <- c("a[1]", "a[2]",
                        "b_trt[1]", "b_trt[2]", "b_trt[3]",
                        "b_prev[1]", "b_prev[2]", "b_prev[3]",
                        "b_time_1", "b_time_2",
                        "b_prev_time[1]", "b_prev_time[2]", "b_prev_time[3]",
                        "b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")
  
  l_spec$non_zero_pars <- c("a[1]", "a[2]",
                            "b_trt[2]", "b_trt[3]",
                            "b_prev[2]", "b_prev[3]",
                            "b_time_1", "b_time_2",
                            "b_prev_time[2]", "b_prev_time[3]",
                            "b_trt_time[2]", "b_trt_time[3]")
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$nsim)
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$nsim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
  
  # convert to the induced expected difference in sojourn times
  l_spec$delta_lab <- c("delta_def", "delta_dis")
  d_sop <- sim18_sop(1:28, l_spec)
  d_tbl <- dcast(d_sop, day ~ trt, value.var = "none")
  
  l_spec$dur_tru <- c( 
    d_tbl[day > 0, sum(soc)],
    d_tbl[day > 0, sum(def)],
    d_tbl[day > 0, sum(dis)]
  )
  names(l_spec$dur_tru) <- l_spec$trt_lab
  
  l_spec$delta_dur_tru <- c( 
    d_tbl[day > 0, sum(def - soc)],
    d_tbl[day > 0, sum(dis - soc)]
  )
  names(l_spec$delta_dur_tru) <- l_spec$delta_lab
  
  l_spec
}

sim18_calibrate_trt <- function(l_spec){
  
  source("R/libs.R")
  source("R/init.R")
  source("R/util.R")
  
  f_cfgsc <- file.path("./etc/sim18/cfg-sim18-v01.yml")
  l_spec <- config::get(file = f_cfgsc)
  
  l_spec <- update_sim18_cfg(l_spec)
  l_spec$t_0 <- seq_along(1:sum(l_spec$N_pt))
  
  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$N_pt)
  
  # where are we at the moment
  l_spec$b_trt
  days <- 1:max(l_spec$visit_days)
  d_sop <- sim18_sop(days, l_spec)
  dcast(d_sop[day > 0], day ~ trt, value.var = "none" )[day %in% c(1, 4, 7, 14, 21, 28)]
  
  
  l_spec$dec_delta_ni <- 1
  explr_interval <- c(-5, 5)
  message("Traget NI margin  : ", l_spec$dec_delta_ni)
  
  f_obj <- function(b_trt) {
    l_spec$b_trt["def"] <- b_trt
    d_sop <- sim18_sop(days, l_spec)
    d_tbl <- dcast(d_sop[day > 0], day ~ trt, value.var = "none")
    delta <- d_tbl[, sum(def - soc)]
    (delta + l_spec$dec_delta_ni)^2
  }
  
  f_trt <- stats::optimise(f = f_obj, interval = explr_interval)
  # nominated treatment effect
  message("Trt effect to induce NI  : ", f_trt$minimum)
  
  # Recompute as a sanity check
  l_spec$b_trt["def"] <- f_trt$minimum
  d_sop <- sim18_sop(days, l_spec)
  
  d_tbl <- dcast(d_sop[day > 0], day ~ trt, value.var = "none" )
  # Duration in no symptom state by trt
  d_tbl[, .(def = sum(def), dis = sum(dis), soc = sum(soc), delta_def = sum(def - soc))]
  
  
}

sim18_ex_mod <- function(){
  
  source("R/libs.R")
  source("R/init.R")
  source("R/util.R")
  
  # CFG
  f_cfgsc <- file.path("./etc/sim18/cfg-sim18-v01.yml")
  l_spec <- config::get(file = f_cfgsc)
  l_spec <- update_sim18_cfg(l_spec)
  l_spec$N_pt <- c(4000, 100, 100)
  l_spec$t_0 <- seq_along(1:sum(l_spec$N_pt)) 
  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$N_pt)
  message("N: ", paste0(cumsum(l_spec$N_pt), collapse = ", "))
  
  # MODEL
  l_dat <- sim18_cohort(l_spec)
  d_cohort <- l_dat$d_cohort
  d_obs <- l_dat$d_obs
  m_1 <- cmdstanr::cmdstan_model("stan/sim18-v01.stan")
  l_mod <- sim18_stan_data(d_obs, l_spec)
  # f_1_optim <- m_1$optimize(data = l_mod, jacobian = TRUE)
  # f_1 <- m_1$laplace(data = l_mod, mode = f_1_optim, draws = 2000, refresh = 0)
  f_1 <- m_1$sample(
    l_mod,
    iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
    parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain,
    refresh = 100, show_exceptions = F,
    max_treedepth = 11
  )
  m_post <- f_1$draws(variables = l_spec$non_zero_pars, format = "matrix")
  d_mu_v_tru <- data.table(
    par = l_spec$non_zero_pars,
    tru = c(l_spec$alpha, l_spec$b_trt[-1],  l_spec$b_prev[-1], 
            l_spec$b_time_1, l_spec$b_time_2, 
            l_spec$b_prev_time[-1], l_spec$b_trt_time[-1]
    ),
    mu = colMeans(m_post), 
    lo = apply(m_post, 2, function(z){quantile(z, prob = 0.025)}),
    hi = apply(m_post, 2, function(z){quantile(z, prob = 0.975)})
  )
  d_mu_v_tru[, cover := lo < tru & tru < hi]
  kableExtra::kbl(d_mu_v_tru[], format = "simple", digits = 4)
  # 
  # par                   tru        mu        lo        hi  cover 
  # ---------------  --------  --------  --------  --------  ------
  #   a[1]              -0.5000   -0.4373   -0.5358   -0.3394  TRUE  
  # a[2]               1.2000    1.2330    1.1314    1.3316  TRUE  
  # b_trt[2]          -0.0368   -0.0145   -0.0921    0.0584  TRUE  
  # b_trt[3]           0.0000    0.0324   -0.0441    0.1067  TRUE  
  # b_prev[2]          2.2000    2.1777    2.0873    2.2623  TRUE  
  # b_prev[3]          1.0000    1.0020    0.9110    1.0890  TRUE  
  # b_time_1          -0.2300   -0.2148   -0.2309   -0.1979  TRUE  
  # b_time_2           0.0030    0.0025    0.0020    0.0030  TRUE  
  # b_prev_time[2]     0.0500    0.0559    0.0470    0.0657  TRUE  
  # b_prev_time[3]     0.0000   -0.0011   -0.0154    0.0126  TRUE  
  # b_trt_time[2]      0.0200    0.0160    0.0067    0.0267  TRUE  
  # b_trt_time[3]      0.0000   -0.0047   -0.0142    0.0054  TRUE  
  
  
  
  d_l_spec <- sim18_cfg_to_dtb(l_spec)
  kableExtra::kbl(d_l_spec[], format = "simple")
  
  # item             val                                                  
  # ---------------  -----------------------------------------------------
  #   p_init           0, 0.4, 0.6                                          
  # nsim             500                                                  
  # mcmc_type        laplace                                              
  # N_pt             4000, 100, 100                                       
  # pt_per_day       1.2                                                  
  # ramp_up_days     60                                                   
  # followup         720                                                  
  # followup_dec     28                                                   
  # visit_days       0, 1, 2, 3, 4, 5, 6, 7, 14, 21, 28                   
  # alpha            -0.5, 1.2                                            
  # b_trt            0, -0.03682139, 0                                    
  # b_prev           0, 2.2, 1                                            
  # b_time_1         -0.23                                                
  # b_time_2         0.003                                                
  # b_trt_time       0, 0.02, 0                                           
  # b_prev_time      0, 0.05, 0                                           
  # dec_delta_ni     1                                                    
  # dec_thresh_ni    0.985                                                
  # dec_thresh_fut   0.8                                                  
  # dur_tru          19.7915123657467, 18.7914187697217, 19.7915123657467 
  # delta_dur_tru    -1.00009359602499, 0     
  
}


sim18_cfg_to_dtb <- function(l_spec){
  
  d_l_spec <- data.table(
    item = c("p_init", "nsim", "mcmc_type", "N_pt", "pt_per_day", 
             "ramp_up_days", "followup", "followup_dec", "visit_days", 
             "alpha", "b_trt", "b_prev", "b_time_1", "b_time_2", 
             "b_trt_time", "b_prev_time", 
             "dec_delta_ni", "dec_thresh_ni", "dec_thresh_fut", 
             "dur_tru", "delta_dur_tru")
  )
  d_l_spec[item == "p_init", val := paste0(l_spec$p_init, collapse = ", ")]
  d_l_spec[item == "nsim", val := paste0(l_spec$nsim, collapse = ", ")]
  d_l_spec[item == "mcmc_type", val := paste0(l_spec$mcmc_type, collapse = ", ")]
  d_l_spec[item == "N_pt", val := paste0(l_spec$N_pt, collapse = ", ")]
  d_l_spec[item == "pt_per_day", val := paste0(l_spec$pt_per_day, collapse = ", ")]
  d_l_spec[item == "ramp_up_days", val := paste0(l_spec$ramp_up_days, collapse = ", ")]
  d_l_spec[item == "followup", val := paste0(l_spec$followup, collapse = ", ")]
  d_l_spec[item == "followup_dec", val := paste0(l_spec$followup_dec, collapse = ", ")]
  d_l_spec[item == "visit_days", val := paste0(l_spec$visit_days, collapse = ", ")]
  d_l_spec[item == "alpha", val := paste0(l_spec$alpha, collapse = ", ")]
  d_l_spec[item == "b_trt", val := paste0(l_spec$b_trt, collapse = ", ")]
  d_l_spec[item == "b_prev", val := paste0(l_spec$b_prev, collapse = ", ")]
  d_l_spec[item == "b_time_1", val := paste0(l_spec$b_time_1, collapse = ", ")]
  d_l_spec[item == "b_time_2", val := paste0(l_spec$b_time_2, collapse = ", ")]
  d_l_spec[item == "b_trt_time", val := paste0(l_spec$b_trt_time, collapse = ", ")]
  d_l_spec[item == "b_prev_time", val := paste0(l_spec$b_prev_time, collapse = ", ")]
  
  d_l_spec[item == "dec_delta_ni", val := paste0(l_spec$dec_delta_ni, collapse = ", ")]
  d_l_spec[item == "dec_thresh_ni", val := paste0(l_spec$dec_thresh_ni, collapse = ", ")]
  d_l_spec[item == "dec_thresh_fut", val := paste0(l_spec$dec_thresh_fut, collapse = ", ")]
  d_l_spec[item == "dur_tru", val := paste0(l_spec$dur_tru, collapse = ", ")]
  d_l_spec[item == "delta_dur_tru", val := paste0(l_spec$delta_dur_tru, collapse = ", ")]
  d_l_spec[]
  
  
  
  
  
}

#' Inhomogeneous Poisson process
#' 
#' Simulate fixed n by thinning
sim_ipp_thinning <- function(
    N = 2500, lambda = 1.52,
    rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}



run_trial <- function(
    ix,
    l_spec,
    return_posterior = F
){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # Get enrolment times for arbitrary large set of patients
  # Simpler to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  
  # events per day
  lambda = l_spec$pt_per_day
  # ramp up over x months 
  rho = function(t) pmin(t/l_spec$ramp_up_days, 1)
  
  # day of enrolment
  l_spec$t_0 <- sim_ipp_thinning(sum(l_spec$N_pt), lambda, rho)
  
  # t_600 <- numeric(100)
  # for(i in 1:100){
  #   t_600[i] = sim_ipp_thinning(sum(l_spec$N_pt), lambda, rho)[600]
  # }
  # hist(t_600)
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ic <- 1 # interim number
  N_analys <- length(l_spec$N_pt)
  
  # summary for model parameters
  d_post_smry_1 <- CJ(
    ic = 1:N_analys,
    par = l_spec$full_pars
  )
  
  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, lo := NA_real_]
  d_post_smry_1[, hi := NA_real_]
  
  # summary sop
  d_post_smry_2 <- CJ(
    ic = 1:N_analys,
    state = l_spec$state_lab,
    trt = l_spec$trt_lab,
    day = 1:max(l_spec$visit_days)
  )
  d_post_smry_2[, mu := NA_real_]
  d_post_smry_2[, lo := NA_real_]
  d_post_smry_2[, hi := NA_real_]
  
  # summary time in state, difference in time in state
  d_post_smry_3 <- CJ(
    ic = 1:N_analys,
    state = l_spec$state_lab,
    par = c(l_spec$trt_lab, l_spec$delta_lab)
  )
  d_post_smry_3[, mu := NA_real_]
  d_post_smry_3[, lo := NA_real_]
  d_post_smry_3[, hi := NA_real_]
  
  # decisions 
  d_post_smry_4 <- CJ(
    ic = 1:N_analys,
    rule = c("ni", "fut"),
    par = l_spec$delta_lab,
    p = NA_real_,
    dec = NA_integer_
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  # full posterior 
  d_post_par <- data.table()
  d_post_sop <- data.table()
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " cohort ", l_spec$ic)
    
    # next chunk of data on pts.
    if(l_spec$ic == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N_pt[l_spec$ic] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N_pt[l_spec$ic] - 1
    }
    
    # id and time
    
    
    # We are assuming that the analysis takes place on pt having reached endpoint
    
    d <- sim18_cohort(l_spec)$d_obs
    
    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # are we at the final analysis or interim?
    if(l_spec$ie == sum(l_spec$N_pt)){
      
      log_info("Trial ", ix, " final analysis, using all pt")
      # only those with events contribute
      l_mod <- sim18_stan_data(dd = copy(d_all), l_spec)
      
    } else {
      # t_0 is the entry time (note that this is repeated for each id if they 
      # have more than one event over the 12 months)
      t_now <- d_all[, max(t_0)]
      
      # we include participants that have had an exacerbation and for whom
      # we have completed the subsequent 90 day follow up
      incl_ids <- d_all[t_0 + l_spec$followup_dec < t_now, unique(id)]
      dd <- copy(d_all[id %in% incl_ids])
      l_mod <- sim18_stan_data(dd, l_spec)
      
    }
    
    stopifnot(l_spec$mcmc_type %in% c("laplace", "hmc"))
    if(l_spec$mcmc_type == "laplace"){
      f_1_optim <- m_1$optimize(data = l_mod, jacobian = TRUE)
      f_1 <- m_1$laplace(data = l_mod, mode = f_1_optim, draws = 2000, refresh = 0)
    } 
    
    if(l_spec$mcmc_type == "hmc"){
      f_1 <- m_1$sample(
        l_mod,
        iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
        parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain,
        refresh = 0, show_exceptions = F,
        max_treedepth = 11
      )
    }
    
    
    # l_mod$y_f <- factor(l_mod$y)
    # f_2 <- MASS::polr(l_mod$y_f ~ l_mod$X, Hess=TRUE)
    # then back transform
    
    m_post <- f_1$draws(variables = l_spec$smry_pars, format = "matrix")
    # colMeans(m_post)
    ii <- 1
    l_spec_mod <- copy(l_spec)
    d_sop <- rbindlist(lapply(1:nrow(m_post), function(ii){
      
      # calculate sop based on parameter estimates at their post mean
      l_spec_mod$alpha <- as.numeric(m_post[ii , c("a[1]", "a[2]")])
      l_spec_mod$b_trt <- as.numeric(m_post[ii , c("b_trt[1]", "b_trt[2]", "b_trt[3]")])  
      l_spec_mod$b_prev <- as.numeric(m_post[ii , c("b_prev[1]", "b_prev[2]", "b_prev[3]")])  
      l_spec_mod$b_time_1 <- as.numeric(m_post[ii , c("b_time_1")])  
      l_spec_mod$b_time_2 <- as.numeric(m_post[ii , c("b_time_2")]) 
      l_spec_mod$b_prev_time <- as.numeric(m_post[ii , c("b_prev_time[1]", "b_prev_time[2]", "b_prev_time[3]")])
      l_spec_mod$b_trt_time <- as.numeric(m_post[ii , c("b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")])  
      
      names(l_spec_mod$b_trt) <- l_spec_mod$trt_lab
      names(l_spec_mod$b_trt_time) <- l_spec_mod$trt_lab
      
      d_tmp <- sim18_sop(days = 1:max(l_spec$visit_days), l_spec_mod)
      d_tmp
      
    }), idcol = "id_draw" )
    
    
    if(return_posterior){
      
      # posterior
      d_post_par <- rbind(
        d_post_par,
        cbind(ic = l_spec$ic, m_post)
      )
      
      # posterior
      d_post_sop <- rbind(
        d_post_sop,
        cbind(ic = l_spec$ic, d_sop)
      )
      
    }
    
    
    d_sop[, gap := day - data.table::shift(day, type = "lag"), keyby = .(id_draw, trt)]
    
    d_sop[, `:=`(
      prop_none = none * gap,
      prop_mild = mild * gap,
      prop_severe = severe * gap)]
    
    
    # days in each state
    d_dur <- melt(
      d_sop[day > 0, .(
        none = sum(none * gap),
        mild = sum(mild * gap),
        severe = sum(severe * gap)
      ), keyby = .(trt, id_draw)], id.vars = c("trt", "id_draw"), 
      variable.name = "state", value.name = "days")
    
    d_dur <- dcast(d_dur, state + id_draw ~ trt, value.var = "days")
    d_dur[, delta_def := def - soc]
    d_dur[, delta_dis := dis - soc]
    
    # sd(d_dur$delta_def)
    # sd(d_dur$delta_dis)
    
    # posterior summary (did the model recover the data generating pars)
    d_post <- data.table(m_post)
    d_post_smry_1[
      data.table(
        ic = l_spec$ic, 
        par = names(d_post),
        mu = colMeans(d_post),
        lo = apply(d_post, 2, function(z){
          quantile(z, prob = 0.025)
        }),
        hi = apply(d_post, 2, function(z){
          quantile(z, prob = 0.975)
        })
      ),
      on = .(ic, par), `:=`(
        mu = i.mu, lo = i.lo, hi = i.hi
      )
    ]
    
    
    # posterior summary - what are the estimated sops
    d_post_smry_2[
      melt(d_sop[day > 0, .(id_draw, trt, day, none, mild, severe)], 
           id.vars = c("id_draw", "trt", "day"), variable.name = "state")[
             , .(ic = l_spec$ic,
                 mu = mean(value), 
                 lo = quantile(value, prob = 0.025),
                 hi = quantile(value, prob = 0.975)),
             keyby = .(state, trt, day)
           ],
      on = .(ic, state, trt, day), `:=`(
        mu = i.mu, lo = i.lo, hi = i.hi
      )
    ]
    
    # posterior summary - days in state and treatment effects
    d_post_smry_3[
      melt(d_dur, id.vars = c("state", "id_draw"), variable.name = "par")[
        , .(
          ic = l_spec$ic,
          mu = mean(value),
          lo = quantile(value, prob  = 0.025),
          hi = quantile(value, prob  = 0.975)
        ), keyby = .(state, par)
      ],
      on = .(ic, par, state), `:=`(
        mu = i.mu, lo = i.lo, hi = i.hi
      )
    ]
    
    
    # posterior summary - evaluating decision rules on the differences
    d_post_smry_4[
      melt(d_dur[
        state == "none", .(id_draw, delta_def, delta_dis)], 
        id.vars = c("id_draw"), variable.name = "par")[
          , .(
            ic = l_spec$ic,
            rule = "ni",
            # probability that difference in number of symptom free days is above 
            # some nominal reduction
            p = mean(value > -l_spec$dec_delta_ni) ,
            dec = NA_real_
          ), keyby = .(par)
        ],
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec
      )
    ]
    d_post_smry_4[ic == l_spec$ic & rule == "ni", dec := as.integer(p > l_spec$dec_thresh_ni)]
    
    d_post_smry_4[
      melt(d_dur[state == "none", .(id_draw, delta_def, delta_dis)], id.vars = c("id_draw"), variable.name = "par")[
        , .(
          ic = l_spec$ic,
          rule = "fut",
          # probability that difference in number of symptom free days is above 
          # some nominal reduction
          p = mean(value < -l_spec$dec_delta_ni) ,
          dec = NA_real_
        ), keyby = .(par)
      ],
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec
      )
    ]
    d_post_smry_4[ic == l_spec$ic & rule == "fut", dec := as.integer(p > l_spec$dec_thresh_fut)]
    
    
    
    # evaluate status - have we stopped?
    d_stop <- d_post_smry_4[
      ic <= l_spec$ic,
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]
    
    if(all(d_stop$resolved)){
      log_info("Stop trial as all questions addressed ", ix)
      stop_enrol <- T
    } else if(any(d_stop$resolved)){
      
      if(d_stop[par == "delta_def", resolved]) {
        l_spec$trt_active["def"] <- FALSE
        
        log_info("Trial ", ix, 
                 " def stopped no further enrolment into this arm, analy id ", 
                 l_spec$ic)
      }
      
      if(d_stop[par == "delta_dis", resolved]) {
        l_spec$trt_active["dis"] <- FALSE
        
        log_info("Trial ", ix, 
                 " dis stopped no further enrolment into this arm, analy id ", 
                 l_spec$ic)
      }
    }
    
    log_info("Trial ", ix, " allocation after cohort ", 
             l_spec$ic, " alloc: ", paste0(l_spec$trt_active, collapse = ", "))
    
    
    # next interim
    l_spec$ic <- l_spec$ic + 1
    
    if(l_spec$ic > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- l_spec$ic - 1
  
  # lobstr::obj_size(d_w)
  # lobstr::obj_size(d_all)
  
  l_ret <- list(
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    d_post_smry_3 = d_post_smry_3,
    d_post_smry_4 = d_post_smry_4,
    d_post_par = d_post_par,
    d_post_sop = d_post_sop,
    stop_at = stop_at,
    l_spec = l_spec
  )
  
  return(l_ret)
}





#' Entry point for running trial simulation in parallel.
run_sim18 <- function(){
  
  tic()
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 3
  }
  
  l_spec <- update_sim18_cfg(l_spec)
  
  return_posterior <- F
  str(l_spec)
  e = NULL
  ix <- 1
  
  ## LOOP -------
  
  
  log_info("Start simulation")
  
  r <- pbapply::pblapply(
    X=1:l_spec$nsim, cl = l_spec$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);
      
      if(ix %in% l_spec$ex_trial_ix){
        return_posterior = T  
      } else {
        return_posterior = F
      }
      
      ll <- tryCatch({
        run_trial(ix, l_spec, return_posterior = return_posterior )
      },
      error=function(e) {
        log_info("ERROR in MCLAPPLY LOOP (see terminal output):")
        message(" ERROR in MCLAPPLY LOOP " , e); message(traceback())
        stop(paste0("Stopping with error ", e))
      })
      
      ll
    })
  
  
  Sys.sleep(10)
  
  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_all
  } ), idcol = "sim")
  
  d_post_smry_1 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_1
  } ), idcol = "sim")
  
  d_post_smry_2 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_2
  } ), idcol = "sim")
  
  d_post_smry_3 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_3
  } ), idcol = "sim")
  
  d_post_smry_4 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_4
  } ), idcol = "sim")
  
  d_post_par <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_par
  } ), idcol = "sim")
  
  d_post_sop <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_sop
  } ), idcol = "sim")
  
  l_spec$toc <- toc()
  l <- list(
    l_spec = l_spec,
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    d_post_smry_3 = d_post_smry_3,
    d_post_smry_4 = d_post_smry_4,
    d_post_par = d_post_par,
    d_post_sop = d_post_sop
  )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim18/sim18-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  message("fname is ", fname)
  qs::qsave(l, file = fname)
  
  
  message("saved")
}






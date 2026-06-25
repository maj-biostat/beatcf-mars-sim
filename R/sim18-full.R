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

l_spec <- list(
  p_init = c(0, 0.4, 0.6),
  desc = "Treatments equivalent",
  nsim = 500,
  mc_cores = 40,
  mcmc_warmup = 1000,
  mcmc_iter = 1000,
  mcmc_chain = 1,
  mcmc_B = 1000,
  nex = 5,
  N_pt  = c( 400, 100, 100),
  pt_per_day = 1.2,
  ramp_up_days = 60,
  followup = 720,
  followup_dec = 28,
  visit_days = c(0, 1, 2, 3, 4, 5, 6, 7, 14, 21, 28),
  state_lab = c("none", "mild", "severe"),
  trt_active = c(1, 1, 1),
  trt_lab = c("soc", "def", "dis"),
  alpha = c(-0.5, 1.2),
  # 0.1418942 -0.4000009 -0.4000000 corresponds to an NI of 0.40
  # 0.1592559 -0.4499971 -0.4500000 corresponds to an NI of 0.45
  # 0.1765427 -0.5000003 -0.5000000 corresponds to an NI of 0.50
  # 0.1937544 -0.5498294 -0.5500000 corresponds to an NI of 0.55 
  b_trt = c(0.0, 0.1592559, 0.0)   ,
  b_prev = c(0, 2, 1),
  b_time_1 = -0.32,
  b_time_2 = 0.003,
  b_gap = c(0, 0),
  b_trt_gap = c(0.0, 0.0, 0.0),
  b_trt_time = c(0, 0, 0),
  dec_delta_ni = 0.5,
  dec_thresh_ni = 0.975,
  dec_thresh_fut = 0.8
  
)


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
  int  ix_gap ;
  int  ix_trt_gap_2;
  int  ix_trt_gap_3;
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
  b ~ normal(0, 1);

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
  b_prev[2] = b[ix_prev_2];
  b_prev[3] = b[ix_prev_3];
  
  real b_time_1 = (b[ix_time_1] / sd_days) - (2.0 * mu_days * b[ix_time_2]) / pow(sd_days, 2);
  real b_time_2 = b[ix_time_2] / pow(sd_days, 2);
  
  vector[2] b_gap;
  b_gap[1] = 0.0;
  b_gap[2] = b[ix_gap];
  
  vector[3] b_trt_gap;
  b_trt_gap[1] = 0.0;
  b_trt_gap[2] = b[ix_trt_gap_2];
  b_trt_gap[3] = b[ix_trt_gap_3];
  
  vector[3] b_trt_time;
  b_trt_time[1] = 0.0;
  b_trt_time[2] = b[ix_trt_time_2] / sd_days;
  b_trt_time[3] = b[ix_trt_time_3] / sd_days;
  
}
"

stan_file <- write_stan_file(stan_model_code)
m_1 <- cmdstanr::cmdstan_model(stan_file)

output_dir_mcmc <- paste0(getwd(), "/tmp")

#' Wrapper to invoke state transition simulation for individual patients used
#' to create cohort
sim18_cohort <- function(l_spec){
  
  id_cohort <- l_spec$is:l_spec$ie
  N_cohort <- length(id_cohort)
  
  # produce data for every day and then chop it down to what we observe
  d_cohort <- data.table(
    id  = l_spec$is:l_spec$ie,
    day = rep(0:l_spec$max_day, N_cohort),
    trt = sample(l_spec$trt_lab[l_spec$trt_active], N_cohort, replace = T)
  )
  setorder(d_cohort, id, day)
  d_cohort[, t_0 := l_spec$t_0[id]]
  
  d_cohort[, state := NA_integer_]
  d_cohort[day==0, state:= sample(l_spec$state_opts, .N, replace = TRUE, prob = l_spec$p_init)]
  
  for(i in seq_len(N_cohort))
  {
    
    rows <- which(d_cohort$id == id_cohort[i])
    
    for(j in 2:length(rows))
    {
      
      prev <- d_cohort$state[rows[j-1]]
      
      # all intervals are structurally 1 because we simulate daily and then 
      # subset to the survey days
      interval <- 1
      
      tt <- d_cohort$day[rows[j]]
      trt <- d_cohort$trt[rows[j]]
      
      lp <-
        # treatment
        l_spec$b_trt[trt] +
        # previous state
        l_spec$b_prev[prev] +
        # time (quadratic)
        l_spec$b_time_1*tt +
        l_spec$b_time_2*tt^2 +
        # time by treatment (linear)
        l_spec$b_trt_time[trt]*tt +
        # gap length - structurally zero and is only relevant in model
        # where we have gaps in observation and it constitutes a nuissance param
        l_spec$b_gap[interval] +
        # trt x gap is similarly structurally zero and is only relevant in 
        # the model to account for differential progression across trt arms
        as.numeric(interval != 1) * l_spec$b_trt_gap[trt] 
      
      
      # not very efficient but it will do
      d_cohort$state[rows[j]] <- sim18_rord_pom_3(lp,l_spec$alpha)
      
    }
    
  }
  
  d_obs <- d_cohort[day %in% l_spec$visit_days]
  
  list(
    d_cohort = d_cohort,
    d_obs = d_obs
  )
  
}



#' Convert sample data.table into lists suitable for stan models
#' 
sim18_stan_data <- function(dd, l_spec){
  
  
  setorder(dd, id, day)
  
  dd[, `:=`(
    prev_state = data.table::shift(state, 1L),
    prev_day   = data.table::shift(day, 1L)
  ), by = id]
  
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
  dd[, x_gap_time := factor(fifelse(dd$gap_len == 1, 1, 2))]
  dd[, x_trt := factor(trt_idx)]
  dd[, x_prev := factor(prev_state, levels = l_spec$state_opts)]
  
  X <- model.matrix(~ x_trt + x_prev +
                      x_time + I(x_time^2) +
                      x_gap_time +
                      x_trt * x_gap_time + 
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
    ix_gap = 7,
    ix_trt_gap_2 = 8,
    ix_trt_gap_3 = 9,
    ix_trt_time_2 = 10,
    ix_trt_time_3 = 11,
    mu_days = mean(dd$day),
    sd_days = sd(dd$day)
  )
  stopifnot(names(X_mod)[ld$ix_trt_2] == "x_trt2")
  stopifnot(names(X_mod)[ld$ix_trt_3] == "x_trt3")
  stopifnot(names(X_mod)[ld$ix_prev_2] == "x_prev2")
  stopifnot(names(X_mod)[ld$ix_prev_3] == "x_prev3")
  stopifnot(names(X_mod)[ld$ix_time_1] == "x_time")
  stopifnot(names(X_mod)[ld$ix_time_2] == "I(x_time^2)")
  stopifnot(names(X_mod)[ld$ix_gap] == "x_gap_time2")
  stopifnot(names(X_mod)[ld$ix_trt_gap_2] == "x_trt2:x_gap_time2")
  stopifnot(names(X_mod)[ld$ix_trt_gap_3] == "x_trt3:x_gap_time2")
  stopifnot(names(X_mod)[ld$ix_trt_time_2] == "x_trt2:x_time")
  stopifnot(names(X_mod)[ld$ix_trt_time_3] == "x_trt3:x_time")
  
  ld
}



# construct 3 state transition probs based on linear predictor and 
# cutpoints from ordinal model
sim18_rord_pom_3 <- function(lp, alpha)
{
  p0 <- plogis(alpha[1] - lp)
  p1 <- plogis(alpha[2] - lp) - p0
  p2 <- 1 - p0 - p1
  
  sample(1:3, 1, prob = c(p0, p1, p2)) 
}

sim18_transition_matrix <- function(day, gap_ix = 1, trt, l_spec)
{
  
  P <- matrix(0,3,3)
  
  for(prev in l_spec$state_opts){
    
    lp <-
      l_spec$b_trt[trt] +
      l_spec$b_prev[prev] +
      l_spec$b_time_1 * day +
      l_spec$b_time_2 * day^2 +
      l_spec$b_trt_time[trt] * day +
      l_spec$b_gap[gap_ix] +
      as.numeric(gap_ix != 1) * l_spec$b_trt_gap[trt]
    
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
      
      if(day == 1){
        gap_ix = 1;
      } 
      
      if(day > 1){
        if(days[i] - days[i-1] == 1){
          gap_ix = 1
        } else {
          gap_ix = 2
        }
      }
      
      Pt <- sim18_transition_matrix(day, gap_ix, trt, l_spec)
      
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
  
  l_spec$b_gap <- unlist(l_spec$b_gap)
  
  l_spec$b_trt_gap <- unlist(l_spec$b_trt_gap)
  names(l_spec$b_trt_gap) <- l_spec$trt_lab
  
  l_spec$p_init <- unlist(l_spec$p_init)
  
  l_spec$smry_pars <- c(
    "a", "b_trt", "b_prev", "b_time_1", "b_time_2", "b_gap", "b_trt_gap", "b_trt_time")
  
  l_spec$full_pars <- c("a[1]", "a[2]",
                        "b_trt[1]", "b_trt[2]", "b_trt[3]",
                        "b_prev[1]", "b_prev[2]", "b_prev[3]",
                        "b_time_1", "b_time_2",
                        "b_gap[1]", "b_gap[2]",
                        "b_trt_gap[1]", "b_trt_gap[2]", "b_trt_gap[3]",
                        "b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")
  
  l_spec$non_zero_pars <- c("a[1]", "a[2]",
                            "b_trt[2]", "b_trt[3]",
                            "b_prev[2]", "b_prev[3]",
                            "b_time_1", "b_time_2",
                            "b_gap[2]",
                            "b_trt_gap[2]", "b_trt_gap[3]",
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
      l_mod <- sim18_stan_data(dd = copy(d_all[id %in% incl_ids]), l_spec)
      
    }
    
    f_1_optim <- m_1$optimize(data = l_mod, jacobian = TRUE)
    f_1 <- m_1$laplace(data = l_mod, mode = f_1_optim, draws = 2000, refresh = 0)
    
    # f_1 <- m_1$sample(
    #   l_mod,
    #   iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
    #   parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain,
    #   refresh = 0, show_exceptions = T,
    #   max_treedepth = 11
    # )
    
    m_post <- f_1$draws(variables = l_spec$smry_pars, format = "matrix")
    ii <- 1
    l_spec_mod <- copy(l_spec)
    d_sop <- rbindlist(lapply(1:nrow(m_post), function(ii){
      
      # calculate sop based on parameter estimates at their post mean
      l_spec_mod$alpha <- as.numeric(m_post[ii , c("a[1]", "a[2]")])
      l_spec_mod$b_trt <- as.numeric(m_post[ii , c("b_trt[1]", "b_trt[2]", "b_trt[3]")])  
      l_spec_mod$b_prev <- as.numeric(m_post[ii , c("b_prev[1]", "b_prev[2]", "b_prev[3]")])  
      l_spec_mod$b_time_1 <- as.numeric(m_post[ii , c("b_time_1")])  
      l_spec_mod$b_time_2 <- as.numeric(m_post[ii , c("b_time_2")]) 
      l_spec_mod$b_gap <- as.numeric(m_post[ii , c("b_gap[1]", "b_gap[2]")])
      l_spec_mod$b_trt_gap <- as.numeric(m_post[ii , c("b_trt_gap[1]", "b_trt_gap[2]", "b_trt_gap[3]")])  
      l_spec_mod$b_trt_time <- as.numeric(m_post[ii , c("b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")])  
      
      names(l_spec_mod$b_trt) <- l_spec_mod$trt_lab
      names(l_spec_mod$b_trt_gap) <- l_spec_mod$trt_lab
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
  
  
  
  # 
  
  return(l_ret)
}





#' Entry point for running trial simulation in parallel.
run_sim18 <- function(){
  
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
  
  
  fname <- paste0("sim18-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  message("fname is ", fname)
  qs::qsave(l, file = fname)
  
  
  message("saved")
}


# Simulation config from first scenario
# > l_spec_1$p_init
# [1] 0.0 0.4 0.6
# > l_spec_1$alpha
# [1] -0.5  1.2
# > l_spec_1$b_trt
# soc       def       dis 
# 0.0000000 0.1592559 0.0000000 
# > l_spec_1$b_prev
# [1] 0 2 1
# > l_spec_1$b_time_1
# [1] -0.32
# > l_spec_1$b_time_2
# [1] 0.003
# > l_spec_1$b_gap
# [1] 0 0
# > l_spec_1$b_trt_time
# soc def dis 
# 0   0   0 
# > l_spec_1$dec_delta_ni
# [1] 0.45
# > l_spec_1$dec_thresh_ni
# [1] 0.975
# > l_spec_1$dec_thresh_fut
# [1] 0.8
# > l_spec_1$visit_days
# [1]  0  1  2  3  4  5  6  7 14 21 28


# Results from simulation - model parameters (posterior mean and 95CrI) by interim
# Source for Table 6.3
# > d_par
# scenario                 desc    ic           par           mu            lo           hi     N
# <int>               <char> <int>        <char>        <num>         <num>        <num> <int>
#   1:        1 Defer at NI Boundary     1          a[1] -0.530256360 -0.8357155031 -0.201161764   400
# 2:        1 Defer at NI Boundary     1          a[2]  1.157840107  0.8471641685  1.482078981   400
# 3:        1 Defer at NI Boundary     1      b_gap[2] -0.783922953 -1.4872662533 -0.148253028   400
# 4:        1 Defer at NI Boundary     1     b_prev[2]  1.896291740  1.7077374534  2.073741823   400
# 5:        1 Defer at NI Boundary     1     b_prev[3]  0.931975811  0.7203691493  1.169913654   400
# 6:        1 Defer at NI Boundary     1      b_time_1 -0.317662611 -0.3803827197 -0.242933414   400
# 7:        1 Defer at NI Boundary     1      b_time_2  0.004160289 -0.0004079629  0.007458714   400
# 8:        1 Defer at NI Boundary     1      b_trt[2]  0.202193574 -0.0915801435  0.517366220   400
# 9:        1 Defer at NI Boundary     1      b_trt[3]  0.053855739 -0.2377208945  0.361441883   400
# 10:        1 Defer at NI Boundary     1 b_trt_time[2] -0.013429627 -0.0625289386  0.035204460   400
# 11:        1 Defer at NI Boundary     1 b_trt_time[3] -0.013982053 -0.0682846255  0.041875638   400
# 12:        1 Defer at NI Boundary     2          a[1] -0.532525721 -0.7981772868 -0.231743873   500
# 13:        1 Defer at NI Boundary     2          a[2]  1.156475320  0.8987212610  1.478950630   500
# 14:        1 Defer at NI Boundary     2      b_gap[2] -0.780046820 -1.4159359865 -0.194380010   500
# 15:        1 Defer at NI Boundary     2     b_prev[2]  1.905003710  1.7339300057  2.068426444   500
# 16:        1 Defer at NI Boundary     2     b_prev[3]  0.938714008  0.7377782631  1.130943746   500
# 17:        1 Defer at NI Boundary     2      b_time_1 -0.320656064 -0.3733071623 -0.258591030   500
# 18:        1 Defer at NI Boundary     2      b_time_2  0.004210258  0.0003516069  0.007169683   500
# 19:        1 Defer at NI Boundary     2      b_trt[2]  0.196100924 -0.0686977579  0.482165297   500
# 20:        1 Defer at NI Boundary     2      b_trt[3]  0.046968943 -0.2232444209  0.308748892   500
# 21:        1 Defer at NI Boundary     2 b_trt_time[2] -0.010954294 -0.0555111007  0.034300791   500
# 22:        1 Defer at NI Boundary     2 b_trt_time[3] -0.011991109 -0.0639537570  0.039633725   500
# 23:        1 Defer at NI Boundary     3          a[1] -0.541418422 -0.7763792703 -0.281673799   600
# 24:        1 Defer at NI Boundary     3          a[2]  1.147886536  0.9032569378  1.431352404   600
# 25:        1 Defer at NI Boundary     3      b_gap[2] -0.779948102 -1.3486130081 -0.227356690   600
# 26:        1 Defer at NI Boundary     3     b_prev[2]  1.909876294  1.7497859112  2.058336277   600
# 27:        1 Defer at NI Boundary     3     b_prev[3]  0.943082967  0.7768105355  1.111060607   600
# 28:        1 Defer at NI Boundary     3      b_time_1 -0.324143531 -0.3713567980 -0.275156704   600
# 29:        1 Defer at NI Boundary     3      b_time_2  0.004324649  0.0009471881  0.006908316   600
# 30:        1 Defer at NI Boundary     3      b_trt[2]  0.190596257 -0.0354300326  0.424845447   600
# 31:        1 Defer at NI Boundary     3      b_trt[3]  0.037112990 -0.2063398726  0.285039400   600
# 32:        1 Defer at NI Boundary     3 b_trt_time[2] -0.009509152 -0.0528950407  0.033302119   600
# 33:        1 Defer at NI Boundary     3 b_trt_time[3] -0.010498981 -0.0651219086  0.035871594   600
# scenario                 desc    ic           par           mu            lo           hi     N
# <int>               <char> <int>        <char>        <num>         <num>        <num> <int>


# Results from simulation - model parameters (posterior mean and 95CrI) by interim
# Posterior parameter summary (derived parameters - sojourn times)
# Source for Table 6.4
# > d_par
# scenario                 desc    ic    par       mu       lo       hi     N
# <int>               <char> <int> <char>    <num>    <num>    <num> <int>
#   1:        1 Defer at NI Boundary     1    def 23.39856 23.03698 23.73489   400
# 2:        1 Defer at NI Boundary     1    dis 23.72880 23.38601 24.05368   400
# 3:        1 Defer at NI Boundary     1    soc 23.66356 23.29854 23.99548   400
# 4:        1 Defer at NI Boundary     2    def 23.40376 23.08081 23.71041   500
# 5:        1 Defer at NI Boundary     2    dis 23.74278 23.39661 24.06334   500
# 6:        1 Defer at NI Boundary     2    soc 23.69119 23.35964 23.98518   500
# 7:        1 Defer at NI Boundary     3    def 23.40948 23.08741 23.71249   600
# 8:        1 Defer at NI Boundary     3    dis 23.75483 23.41878 24.05641   600
# 9:        1 Defer at NI Boundary     3    soc 23.70399 23.39250 23.96202   600

# Results from simulation - model parameters (posterior mean and 95CrI) by interim
# Posterior parameter summary (derived parameters - difference in sojourn times)
# Source for Table 6.5
# > d_par
# scenario                 desc    ic       par          mu         lo        hi     N
# <int>               <char> <int>    <char>       <num>      <num>     <num> <int>
#   1:        1 Defer at NI Boundary     1 delta_def -0.26500580 -0.6834766 0.1765732   400
# 2:        1 Defer at NI Boundary     1 delta_dis  0.06523406 -0.4120628 0.4977506   400
# 3:        1 Defer at NI Boundary     2 delta_def -0.28742956 -0.6985776 0.1765732   500
# 4:        1 Defer at NI Boundary     2 delta_dis  0.05159393 -0.3415595 0.4666675   500
# 5:        1 Defer at NI Boundary     3 delta_def -0.29450587 -0.7025035 0.1765732   600
# 6:        1 Defer at NI Boundary     3 delta_dis  0.05083606 -0.3080355 0.4120981   600



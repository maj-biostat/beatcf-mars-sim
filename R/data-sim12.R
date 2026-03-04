

# source dependencies
toks <- unlist(data.table::tstrsplit(getwd(), "/")) 
if(toks[length(toks)] == "beatcf-mars-sim"){
  prefix_cfg <- "./etc/sim12/"
  prefix_stan <- "./stan"
  prefix_fig <- "./fig"
  prefix_data <- "./data"
  prefix_r <- "./R"
} else {
  prefix_cfg <- "../etc/sim12/"
  prefix_stan <- "../stan"
  prefix_fig <- "../fig"
  prefix_data <- "../data"
  prefix_r <- "../R"
}

source(paste0(prefix_r, '/libs.R'))
source(paste0(prefix_r, '/init.R'))
source(paste0(prefix_r, '/util.R'))



# todo 
# review betancourt warping paper etc


#' State transitions for patient to followup time
#' 
get_sim12_pt <- function(
    l_spec
  ){
  
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    l_spec$age_min, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd, lower.tail = T)
  p_gt_age_upr <- plnorm(
    l_spec$age_max, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(1, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  age <- qlnorm(u, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd)
  
  
  # standardised baseline ppfev and convert so that unit change actually 
  # corresponds to a 10% shift
  ppfev_baseline <- (ppfev_0(age, sd_ppfev = 3) - l_spec$ppfev_ref) / l_spec$ppfev_increment
  
  # indep gamma frailty but with same param values for each transition
  u_he <- rgamma(1, shape = l_spec$g_a, rate = l_spec$g_r)
  u_eh <- rgamma(1, shape = l_spec$g_a, rate = l_spec$g_r)
  
  # entry time is zero, in a healthy state
  t <- 0
  state <- "H"
  states_list <- list()
  exac_count  <- 0
  
  names(l_spec$b_trt) <- l_spec$trt_lab
  
  while (t < l_spec$followup) {
    
    if (state == "H") {
      
      # draw gap time H->E
      lambda <- u_he * exp(l_spec$mu_exacerb + 
                      l_spec$b_ppfev_exacerb * ppfev_baseline)
      
      gap <- flexsurv::rweibullPH(1, shape = l_spec$shape_he, scale = lambda)
      
      t_next <- t + gap
      if (t_next > l_spec$followup) {
        states_list[[length(states_list)+1]] <-
          data.table(state="H", start=t, stop=l_spec$followup, trt="none")
        break
      }
      
      states_list[[length(states_list)+1]] <-
        data.table(state="H", start=t, stop=t_next, trt="none")
      
      state <- "E"
      t <- t_next
      
    } else {
      
      # increment number of exacerbations
      exac_count <- exac_count + 1
      
      # we randomise treatment at the start of the exacerbation we just hit
      trt_options <- l_spec$trt_lab[l_spec$trt_active]
      
      trt <- sample(trt_options,1)
      # and then compute the recovery time on that basis
      lambda <- u_eh * exp(l_spec$mu_recov + 
                      l_spec$b_ppfev_recov * ppfev_baseline + 
                      l_spec$b_trt[trt])
      gap <- flexsurv::rweibullPH(1, shape = l_spec$shape_eh, scale = lambda)
      
      t_next <- t + gap
      
      if (t_next > l_spec$followup) {
        states_list[[length(states_list)+1]] <-
          data.table(state="E", start=t, stop=l_spec$followup, trt=trt)
        break
      }
      
      states_list[[length(states_list)+1]] <-
        data.table(state="E", start=t, stop=t_next, trt=trt)
      
      # back to healthy
      state <- "H"
      t <- t_next
      
    }
    
  }
  
  # the last observation should probably be censored because we do not see
  # the completion of the exacerbation or healthy spell
  d_pt <- rbindlist(states_list)
  d_pt <- cbind(age = age, ppfev_baseline = ppfev_baseline, d_pt)
  d_pt[, age := age + start/l_spec$followup]
  d_pt
  
}

#' Wrapper to invoke state transition simulation for individual patients used
#' to create cohort
#' 
get_sim12_cohort <- function(l_spec){
  
  pt_list <- list()
  
  id_cohort <- l_spec$is:l_spec$ie
  i <- 1
  for(i in seq_along(id_cohort)){
    
    pt_list[[i]] <- cbind(id = id_cohort[i],  t0 = l_spec$t0[i], get_sim12_pt(l_spec))
    
    i <- i + 1
  }
  
  d_cohort <- rbindlist(pt_list)
  d_cohort[, dur := stop - start]
  d_cohort[]
  
}


#' Convert sample data.table into lists suitable for stan models
#' 
get_sim12_stan_data <- function(d_all){
  
  # The analysis proceeds on only those that have exacerbations
  d_tmp <- copy(d_all[state == "E"])
  d_tmp[, ppfev_std := ppfev_baseline]
  d_tmp[, defer := as.numeric(trt == "defer")]
  d_tmp[, discont := as.numeric(trt == "discont")]
  
  # However, a sojourn is censored if the patient was still in this state 
  # at followup end i.e. we never see them recover and as such their stop time 
  # will equal the followup duration exactly
  d_tmp[, censored := as.numeric(stop == l_spec$followup)]
  
  # rejig the id indexs so that they are contiguous
  d_tmp[, id_idx := as.integer(factor(id))]
  N_id <- d_tmp[, uniqueN(id_idx)]
  
  # and now we have two datasets
  # Completed sojourns: contribute pdf
  d_obs  <- d_tmp[censored == 0]
  # Censored sojourns: contribute survival function only
  d_cens <- d_tmp[censored == 1]
  
  # Every subject who appears in d_tmp appears in d_obs, d_cens, or both.
  # d_cens contains subjects whose final spell was a censored E sojourn.
  # d_obs contains subjects who completed at least one E->H transition.
  # Neither set is guaranteed to contain all subjects in d_tmp.
  # Subjects who were only observed in a H state do not appear in either
  # d_obs or d_cens, as they have no contribution to the EH likelihood.
  
  
  X_obs <- model.matrix(~-1 + ppfev_std + defer + discont, data = d_obs)
  X_cens <- model.matrix(~-1 + ppfev_std + defer + discont, data = d_cens)
  
  ld_eh <- list(
    # completed sojourns
    N_obs   = nrow(d_obs),
    y_obs   = d_obs$dur,
    id_obs  = d_obs$id_idx,
    X_obs   = X_obs,
    # censored final sojourns
    N_cens  = nrow(d_cens),
    y_cens  = d_cens$dur,
    id_cens = d_cens$id_idx,
    X_cens  = X_cens,
    # dimensions
    N_id    = N_id,
    P       = ncol(X_obs),
    compute_rmst = 1,
    tau_rmst = l_spec$rmst_eh_horizon,
    trt_defer_col = 2, trt_discont_col = 3,
    pri_s_u = l_spec$pri_s_u
  )
  
  # Every subject appears in d_tmp_he (state == "H") since all subjects
  # start in H and must have at least one H spell.
  d_tmp <- copy(d_all[state == "H"])
  d_tmp[, ppfev_std := ppfev_baseline]
  
  # However, a sojourn is censored if the patient was still in this state 
  # at followup end i.e. we never see them transition to E and as such their 
  # stop time will equal the followup duration exactly
  d_tmp[, censored := as.numeric(stop == l_spec$followup)]
  
  # rejig the id indexs so that they are contiguous - this shouldn't be an issue
  # here but am doing it anyway so that we create the same id_ix variable to 
  # use in the stan model
  d_tmp[, id_idx := as.integer(factor(id))]
  N_id <- d_tmp[, uniqueN(id_idx)]
  
  # and now we have two datasets
  # Completed sojourns: contribute pdf
  d_obs  <- d_tmp[censored == 0]
  # Censored sojourns: contribute survival function only
  d_cens <- d_tmp[censored == 1]
  
  X_obs <- model.matrix(~-1 + ppfev_std, data = d_obs)
  X_cens <- model.matrix(~-1 + ppfev_std, data = d_cens)
  
  ld_he <- list(
    # completed sojourns
    N_obs   = nrow(d_obs),
    y_obs   = d_obs$dur,
    id_obs  = d_obs$id_idx,
    X_obs   = X_obs,
    # censored final sojourns
    N_cens  = nrow(d_cens),
    y_cens  = d_cens$dur,
    id_cens = d_cens$id_idx,
    X_cens  = X_cens,
    # dimensions
    N_id    = N_id,
    P       = ncol(X_obs),
    compute_rmst = 0,
    tau_rmst = l_spec$rmst_he_horizon,
    trt_defer_col = 2, trt_discont_col = 3,
    pri_s_u = l_spec$pri_s_u
  )
  
  list(
    ld_eh = ld_eh,
    ld_he = ld_he
  )
  
}


#' Simulate exacerbations for artificial sample populations of arbitrary size 
#' based on the posterior parameter draws and reference data set for 
#' the sample covariate distribution. Used to calculate probability of 
#' recovery over time and rmst for the exacerbation state.
#' 
#' For each posterior draw, obtain the durations times for exacerbation for a
#' population of patients characterised by the covariate specification and
#' nominated treatment assignment.
#' Compute the proportion recovered over the evaluation interval to the 
#' max followup period stipulated in t_window
#' The proportion recovered is our proxy for the probability of recovered
#' with uncertainty propagated from the posterior draws.
#' Compute the rmst (to t_window) as the average of the exacerbation duration
sim_episode <- function(
    # posterior draws (EH model only needed for primary metric)
  po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
  po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
  # episode window for recovery metrics 
  # (for bound on restricted mean time to recovery)
  t_window   = 30,
  # evaluate at each time increment
  eval_times = seq(0, 30, by = 1),
  trt = "soc",
  # sample data
  d_cohort
  
) {
  
  S <- length(po_eh_shape)
  N_pop <- nrow(d_cohort)
  n_eval  <- length(eval_times)
  n_draws <- length(po_eh_shape)
  
  # storage: recovery probability at each eval time, per draw
  # and RMST per draw
  p_recovered <- matrix(NA_real_, nrow = n_eval, ncol = S)
  rmst        <- numeric(S)
  
  s <- 1
  for (s in seq_len(S)) {
    
    eh_shape <- po_eh_shape[s]
    eh_ua    <- po_eh_u_a[s]
    eh_btrt  <- switch(trt,
                       soc   = 0,
                       defer = po_eh_b_defer[s],
                       discont = po_eh_b_discont[s]
    )
    
    # scale parameter: mu = exp(b0 + b_ppfev * ppfev_std + b_trt)
    log_mu <- po_eh_b_0[s] + po_eh_b_ppfev[s] * d_cohort$ppfev_baseline + eh_btrt
    mu     <- exp(log_mu)
    
    # draw frailty for each patient, then scale
    u_i   <- rgamma(N_pop, shape = eh_ua, rate = eh_ua)
    lam_i <- u_i * mu   # patient-specific Weibull-PH scale
    
    # simulate sojourn time in E for each patient
    # WeibullPH: S(t) = exp(-lambda * t^alpha)
    # Inverse CDF: t = (-log(U) / lambda)^(1/alpha)
    soj <- (-log(runif(N_pop)) / lam_i)^(1 / eh_shape)
    
    # recovery probability at each eval time:
    # proportion who have recovered (soj <= t) by time t
    # vectorised over eval_times
    p_recovered[, s] <- colMeans(outer(soj, eval_times, "<="))
    
    # RMST: expected days recovered-free within window
    # = mean of min(sojourn, t_window)
    rmst[s] <- mean(pmin(soj, t_window))
    
  }
  
  # summarise recovery probability curve
  d_recovery <- data.table(
    time   = eval_times,
    mean   = rowMeans(p_recovered),
    lo     = apply(p_recovered, 1, quantile, 0.025),
    hi     = apply(p_recovered, 1, quantile, 0.975)
  )
  
  # summarise RMST
  d_rmst <- data.table(
    rmst_mu = mean(rmst),
    rmst_lo   = quantile(rmst, 0.025),
    rmst_hi   = quantile(rmst, 0.975)
  )
  
  
  list(recovery_curve = d_recovery, rmst = d_rmst)
}


#' Simulate sojourn times for h-> and e->h for an arbitrary sized population 
#' representative of the sample distribution based on the observed covariate
#' distribution up to the maximum follow up time.
#' 
#' Used to estimate the total number of exacerbations under each treatment 
#' group as well as the the mean total duration in an exacerbations state 
#' over the followup period with uncertainty is propogated from the posterior.
sim_trajectory <- function(
    # HE posterior draws
  po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
  # EH posterior draws
  po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
  po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
  ix_sampl,
  # settings
  followup  = 365,
  trt       = "soc",
  max_trans = 50L,
  d_cohort
) {
  
  S <- length(ix_sampl)
  N_pop <- nrow(d_cohort)
  
  # per-draw summaries
  mean_time_E  <- numeric(S)
  mean_n_exac  <- numeric(S)
  
  s <- 1
  for (s in seq_len(S)) {
    
    # parameters
    he_shape <- po_he_shape[ix_sampl[s]]
    he_ua    <- po_he_u_a[ix_sampl[s]]
    he_mu    <- exp(po_he_b_0[ix_sampl[s]] + po_he_b_ppfev[ix_sampl[s]] * d_cohort$ppfev_baseline)
    
    eh_btrt  <- switch(trt,
                       soc   = 0,
                       defer = po_eh_b_defer[s],
                       discont = po_eh_b_discont[s]
    )
    eh_shape <- po_eh_shape[ix_sampl[s]]
    eh_ua    <- po_eh_u_a[ix_sampl[s]]
    eh_mu    <- exp(po_eh_b_0[ix_sampl[s]] + po_eh_b_ppfev[ix_sampl[s]] * d_cohort$ppfev_baseline + eh_btrt)
    
    # frailties
    u_he <- rgamma(N_pop, he_ua, he_ua)
    u_eh <- rgamma(N_pop, eh_ua, eh_ua)
    
    # patient-specific scales
    lam_he <- u_he * he_mu
    lam_eh <- u_eh * eh_mu
    
    # simulate all sojourn times in bulk
    # rows = patients, cols = transitions
    n_he <- ceiling(max_trans / 2L)
    n_eh <- floor(max_trans / 2L)
    
    # everyone starts in H state.
    # H sojourns (odd transitions: 1, 3, 5, ...) period in healthy state
    # simulating from weibull distribution in a vectorised way by recycling
    # lam_he (it saves me from looping over lam_he)
    # Recycling: say if the first term were 1:100 and the second were 1:10, we
    # would get:
    # 1/1, 2/2, 3/3, ..., 10/10, 
    # 11/1, 12/2, 13/3, ..., 20/10, ...
    # 91/1, 92/2, 93/3, ..., 98/8, 99/9, 100/10
    soj_he <- matrix(
      (-log(runif(N_pop * n_he)) / lam_he)^(1 / he_shape),
      nrow = N_pop, ncol = n_he
    )
    
    # E sojourns (even transitions: 2, 4, 6, ...)
    soj_eh <- matrix(
      (-log(runif(N_pop * n_eh)) / lam_eh)^(1 / eh_shape),
      nrow = N_pop, ncol = n_eh
    )
    
    # combine the durations by interleaving them into an N_pop x max_trans 
    # matrix - note that this will almost certainly contain more transitions 
    # than are relevant to our follow up period
    # each row contains a pt trajectory of sojourn (gap) times
    soj_all <- matrix(NA_real_, nrow = N_pop, ncol = n_he + n_eh)
    soj_all[, seq(1, by = 2, length.out = n_he)] <- soj_he
    soj_all[, seq(2, by = 2, length.out = n_eh)] <- soj_eh
    # soj_all[1:10, 1:10]
    
    # cumulative transition times
    cum_times <- t(apply(soj_all, 1, cumsum))
    # cum_times[1:10, 1:10]
    
    # cap at followup: entry time into each spell
    # entry to spell k = cum_times[, k-1] (entry to spell 1 = 0)
    entry_times <- cbind(0, cum_times[, -ncol(cum_times)])
    exit_times  <- cum_times
    
    # actual time spent in each spell = min(exit, followup) - min(entry, followup)
    # spells starting after followup contribute 0
    time_in_spell <- pmax(
      pmin(exit_times, followup) - pmin(entry_times, followup),
      0
    )
    # sanity, this should be equal to follow up for all rows
    # rowSums(time_in_spell[1:50, ])
    
    # E spells are even-indexed columns (2, 4, 6, ...)
    e_cols <- seq(2, by = 2, length.out = n_eh)
    # total time in E state
    # time_in_spell[, e_cols, drop = FALSE][1:5, 1:10]
    time_E <- rowSums(time_in_spell[, e_cols, drop = FALSE])
    
    # number of exacerbations = number of E spells that actually started
    # i.e. entry time < followup for even-indexed columns
    n_exac <- rowSums(entry_times[, e_cols, drop = FALSE] < followup)
    
    # diagnostic: warn if cap was hit
    if (any(cum_times[, max_trans] < followup))
      warning(sprintf(
        "Draw %d: %d patients hit max_trans cap - increase max_trans",
        s, sum(cum_times[, max_trans] < followup)
      ))
    
    mean_time_E[s] <- mean(time_E)
    mean_n_exac[s] <- mean(n_exac)
    
  }
  
  list(
    time_E = data.table(
      tot_t_mu = mean(mean_time_E),
      tot_t_lo   = quantile(mean_time_E, 0.025),
      tot_t_hi   = quantile(mean_time_E, 0.975)
    ),
    n_exac = data.table(
      n_exac_mu = mean(mean_n_exac),
      n_exac_lo   = quantile(mean_n_exac, 0.025),
      n_exac_hi   = quantile(mean_n_exac, 0.975)
    ),
    # keep draw-level summaries for contrasts
    draws = data.table(
      draw       = seq_len(S),
      mean_time_E = mean_time_E,
      mean_n_exac = mean_n_exac
    )
  )
}



#' Wraps sim_episode to compute probability of recovery and rmst for 
#' exacerbation state under each treatment arm.
#' 
compare_trts <- function(
    # arms to compare
    arms = c("soc", "defer", "discont"), 
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    # episode window for recovery metrics
    t_window   = 30,
    eval_times = seq(0, 30, by = 1),
    # sample data
    d_cohort
    ) {
  
  
  trt <- arms[1]
  results <- lapply(arms, function(trt) {
    
    res <- sim_episode(
      po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
      po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
      t_window, 
      eval_times,
      trt = trt, # simulate each arm independently
      d_cohort   # sample data
      )
    res$recovery_curve[, trt := trt]
    list(curve = res$recovery_curve, rmst = cbind(trt = trt, res$rmst))
  })
  
  curves <- rbindlist(lapply(results, `[[`, "curve"))
  rmsts  <- rbindlist(lapply(results, `[[`, "rmst"))
  
  # pairwise RMST contrasts (defer vs soc, discont vs soc)
  # re-run with paired draws to get proper posterior contrast
  list(curves = curves, rmsts = rmsts)
}







#' Computes rmst across treatment groups based on draws from the posterior 
#' and the covariate distribution using the sample data.
compare_rmst <- function(
    trt_a = "soc", trt_b = "defer",
    
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    # episode window for recovery metrics
    t_window   = 30,
    delta_ni = 0,
    # sample data for covariate profile
    d_cohort
) {
  
  
  S <- length(po_eh_shape)
  N_pop <- nrow(d_cohort)
  
  rmst_a <- rmst_b <- numeric(S)
  
  for (s in seq_len(S)) {
    
    for (trt in c(trt_a, trt_b)) {
      eh_btrt <- switch(trt, soc = 0,
                        defer = po_eh_b_defer[s],
                        discont = po_eh_b_discont[s])
      mu  <- exp(po_eh_b_0[s] + po_eh_b_ppfev[s] * d_cohort$ppfev_baseline + eh_btrt)
      u_i <- rgamma(N_pop, po_eh_u_a[s], po_eh_u_a[s])
      soj <- (-log(runif(N_pop)) / (u_i * mu))^(1 / po_eh_shape[s])
      val <- mean(pmin(soj, t_window))
      if (trt == trt_a) rmst_a[s] <- val else rmst_b[s] <- val
    }
  }
  
  contrast <- rmst_b - rmst_a
  data.table(
    trt_a = trt_a,
    trt_b = trt_b,
    delta_mu = mean(contrast),
    delta_lo = quantile(contrast, 0.025),
    delta_hi = quantile(contrast, 0.975),
    # does trt_b result in a rmst that is above the rmst for trt_a by less
    # than the ni margin
    pr_lt_ni = mean(contrast < delta_ni)   
  )
}

#' Wrapper for simulating trajectory under arbitrary population.
#' 
compare_trajectory <- function(
    trt_a = "soc",
    trt_b = "defer",
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    # indexes the posterior draws to use
    ix_sampl,
    followup  = 365,
    max_trans = 300L,
    d_cohort
) {
  # pass same draw_idx to both by fixing the RNG seed approach:
  # easier to just run both inside one loop with shared draw_idx
  
  res_a <- sim_trajectory(
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    ix_sampl,
    followup, 
    trt = trt_a,
    max_trans,
    d_cohort
  )
  
  res_b <- sim_trajectory(
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    ix_sampl,
    followup, 
    trt = trt_b,
    max_trans,
    d_cohort
  )
  
  contrast_time  <- res_b$draws$mean_time_E - res_a$draws$mean_time_E
  contrast_nexac <- res_b$draws$mean_n_exac - res_a$draws$mean_n_exac
  
  list(
    time_E = data.table(
      trt_a = trt_a,  
      trt_b = trt_b,
      mean_E_trt_a = res_a$time_E$mean,
      mean_E_trt_b = res_b$time_E$mean,
      delta_mu = mean(contrast_time),
      delta_lo = quantile(contrast_time, 0.025),
      delta_hi = quantile(contrast_time, 0.975),
      pr_positive = mean(contrast_time > 0)
    ),
    n_exac = data.table(
      trt_a = trt_a,  trt_b = trt_b,
      n_exac_trt_a = res_a$n_exac$mean,
      n_exac_trt_b = res_b$n_exac$mean,
      delta_mu = mean(contrast_nexac),
      delta_lo = quantile(contrast_nexac, 0.025),
      delta_hi = quantile(contrast_nexac, 0.975),
      pr_positive = mean(contrast_nexac > 0)
    )
  )
}




#' Prototyping.
example_stan <- function(){
  
  
  m3 <- cmdstanr::cmdstan_model("stan/sim12-v03.stan")
  
  # some of the cohort will have no exacerbation and so the event indicator 
  # for the healthy state record should have an event indicator set to zero 
  # due to administrative censoring
  
  # similarly, the last record for each pt will almost certainly be censored
  # again, this is administrative censoring at 365 days
  
  l_spec <- get_demo_spec()
  
  # l_spec$shape_he <- 2.9
  l_spec$mu_exacerb <- -15.6
  d_cohort <- get_sim12_cohort(l_spec)
  
  weibullPH_summary_stats(shape_ph = l_spec$shape_he, scale_ph = exp(l_spec$mu_exacerb))
  # mean
  exp(l_spec$mu_exacerb)^(-1/l_spec$shape_he) * gamma(1 + 1/l_spec$shape_he)
  
  hist(flexsurv::rweibullPH(1e5, shape = l_spec$shape_he, scale = exp(l_spec$mu_exacerb)))
  
  d_tbl <- d_cohort[, .N, by = .(id, state)]
  all_ids <- unique(d_cohort$id); all_states <- c("H", "E")
  
  d_grid <- CJ(id = all_ids, state = all_states)
  d_tbl_full <- d_tbl[d_grid, on = .(id, state)]
  d_tbl_full[is.na(N), N := 0]
  # if a pt has one exacerbation but is censored (ie we never see the 
  # recovery time because it exceeds fu) then the below summary will show
  # a min of 1 for both E and H
  table(d_tbl_full$state, d_tbl_full$N)
  
  
  
  # Exacerbations
  # The analysis proceeds on only those that have exacerbations
  d_tmp <- copy(d_cohort[state == "E"])
  d_tmp[, ppfev_std := ppfev_baseline]
  d_tmp[, defer := as.numeric(trt == "defer")]
  d_tmp[, discont := as.numeric(trt == "discont")]
  
  # However, a sojourn is censored if the patient was still in this state 
  # at followup end i.e. we never see them recover and as such their stop time 
  # will equal the followup duration exactly
  d_tmp[, censored := as.numeric(stop == l_spec$followup)]
  
  # rejig the id indexs so that they are contiguous
  d_tmp[, id_idx := as.integer(factor(id))]
  N_id <- d_tmp[, uniqueN(id_idx)]
  
  # and now we have two datasets
  # Completed sojourns: contribute pdf
  d_obs  <- d_tmp[censored == 0]
  # Censored sojourns: contribute survival function only
  d_cens <- d_tmp[censored == 1]
  
  # Every subject who appears in d_tmp appears in d_obs, d_cens, or both.
  # d_cens contains subjects whose final spell was a censored E sojourn.
  # d_obs contains subjects who completed at least one E->H transition.
  # Neither set is guaranteed to contain all subjects in d_tmp.
  # Subjects who were only observed in a H state do not appear in either
  # d_obs or d_cens, as they have no contribution to the EH likelihood.
  
  
  X_obs <- model.matrix(~-1 + ppfev_std + defer + discont, data = d_obs)
  X_cens <- model.matrix(~-1 + ppfev_std + defer + discont, data = d_cens)
  
  ld_1 <- list(
    # completed sojourns
    N_obs   = nrow(d_obs),
    y_obs   = d_obs$dur,
    id_obs  = d_obs$id_idx,
    X_obs   = X_obs,
    # censored final sojourns
    N_cens  = nrow(d_cens),
    y_cens  = d_cens$dur,
    id_cens = d_cens$id_idx,
    X_cens  = X_cens,
    # dimensions
    N_id    = N_id,
    P       = ncol(X_obs),
    compute_rmst = 1,
    tau_rmst = l_spec$rmst_eh_horizon,
    trt_defer_col = 2, trt_discont_col = 3,
    pri_s_u = l_spec$pri_s_u
  )
  
  # Every subject appears in d_tmp_he (state == "H") since all subjects
  # start in H and must have at least one H spell.
  d_tmp <- copy(d_cohort[state == "H"])
  d_tmp[, ppfev_std := ppfev_baseline]
  
  # However, a sojourn is censored if the patient was still in this state 
  # at followup end i.e. we never see them transition to E and as such their 
  # stop time will equal the followup duration exactly
  d_tmp[, censored := as.numeric(stop == l_spec$followup)]
  
  # rejig the id indexs so that they are contiguous - this shouldn't be an issue
  # here but am doing it anyway so that we create the same id_ix variable to 
  # use in the stan model
  d_tmp[, id_idx := as.integer(factor(id))]
  N_id <- d_tmp[, uniqueN(id_idx)]
  
  # and now we have two datasets
  # Completed sojourns: contribute pdf
  d_obs  <- d_tmp[censored == 0]
  # Censored sojourns: contribute survival function only
  d_cens <- d_tmp[censored == 1]
  
  # Every subject who appears in d_tmp appears in d_obs, d_cens, or both.
  # d_cens contains subjects whose final spell was a censored E sojourn.
  # d_obs contains subjects who completed at least one E->H transition.
  # Neither set is guaranteed to contain all subjects in d_tmp.
  # Subjects who were only observed in a H state do not appear in either
  # d_obs or d_cens, as they have no contribution to the EH likelihood.
  
  X_obs <- model.matrix(~-1 + ppfev_std, data = d_obs)
  X_cens <- model.matrix(~-1 + ppfev_std, data = d_cens)
  
  ld_2 <- list(
    # completed sojourns
    N_obs   = nrow(d_obs),
    y_obs   = d_obs$dur,
    id_obs  = d_obs$id_idx,
    X_obs   = X_obs,
    # censored final sojourns
    N_cens  = nrow(d_cens),
    y_cens  = d_cens$dur,
    id_cens = d_cens$id_idx,
    X_cens  = X_cens,
    # dimensions
    N_id    = N_id,
    P       = ncol(X_obs),
    compute_rmst = 0,
    tau_rmst = 365,
    trt_defer_col = 2, trt_discont_col = 3,
    pri_s_u = l_spec$pri_s_u
  )
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  foutname <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"), "sim12")
  
  # # mcmc
  f_1_1 <- m3$sample(
    ld_1, iter_warmup = 500, iter_sampling = 500,
    parallel_chains = 1, chains = 1, refresh = 10, show_exceptions = T,
    max_treedepth = 11,
    output_dir = output_dir_mcmc,
    output_basename = foutname
  )
  f_1_1$summary(variables = c("shape", "b_0", "b", "u_a", "rmst", "delta"))
  
  f_1_2 <- m3$sample(
    ld_2, iter_warmup = 500, iter_sampling = 500,
    parallel_chains = 1, chains = 1, refresh = 10, show_exceptions = T,
    max_treedepth = 11,
    output_dir = output_dir_mcmc,
    output_basename = foutname
  )
  f_1_2$summary(variables = c("shape", "b_0", "b", "u_a", "rmst", "delta"))
  # c(l_spec$shape_he, l_spec$shape_eh)
  # c(l_spec$mu_exacerb, l_spec$mu_recov)
  # c(l_spec$b_trt)
  # l_spec$b_ppfev_recov
  
  # laplace approx
  f_2_1_optim <- m3$optimize(data = ld_1, jacobian = TRUE)
  f_2_1 <- m3$laplace(data = ld_1, mode = f_2_1_optim, draws = 2000)
  f_2_1$summary(variables = c("shape", "b_0", "b", "u_a", "rmst", "delta"))
  c(l_spec$shape_eh, l_spec$mu_recov, l_spec$b_ppfev_recov, l_spec$b_trt[-1], l_spec$g_a)
  
  f_2_2_optim <- m3$optimize(data = ld_2, jacobian = TRUE)
  f_2_2 <- m3$laplace(data = ld_2, mode = f_2_2_optim, draws = 2000)
  f_2_2$summary(variables = c("shape", "b_0", "b", "u_a"))
  c(l_spec$shape_he, l_spec$mu_exacerb, l_spec$b_ppfev_exacerb, l_spec$g_a)
  
  
  
  # Posterior predictive simulation
  d_post_eh <- data.table(
    f_1_1$draws(
      format = "matrix", 
      variables = c("shape", "b_0", "b", "u_a"))
  )
  
  po_eh_shape <- d_post_eh$shape
  po_eh_b_0 <- d_post_eh$b_0
  po_eh_b_ppfev <- d_post_eh$`b[1]`
  po_eh_b_defer <- d_post_eh$`b[2]`
  po_eh_b_discont <- d_post_eh$`b[3]`
  po_eh_u_a <- d_post_eh$u_a
  
  d_post_he <- data.table(
    f_1_2$draws(
      format = "matrix", 
      variables = c("shape", "b_0", "b", "u_a"))
  )
  
  po_he_shape <- d_post_he$shape
  po_he_b_0 <- d_post_he$b_0
  po_he_b_ppfev <- d_post_he$`b[1]`
  po_he_u_a <- d_post_he$u_a
  
  arms = c("soc", "defer", "discont")
  # episode window for recovery metrics
  t_window   = l_spec$rmst_eh_horizon
  eval_times = seq(0, l_spec$rmst_eh_horizon, by = 1)
  
  
  
  l_res_1 <- compare_trts(
    arms = c("soc", "defer", "discont"),
    po_eh_shape, 
    po_eh_b_0, po_eh_b_ppfev, po_eh_b_defer, po_eh_b_discont, 
    po_eh_u_a,
    # episode window for recovery metrics
    t_window   = l_spec$rmst_eh_horizon,
    eval_times = seq(0, l_spec$rmst_eh_horizon, by = 1),
    # sample data for covariate profile
    d_cohort
  )
  l_res_1
  
  
  p1 <- ggplot(l_res_1$curves, aes(x = time, y = mean, group = trt, col = trt)) +
    geom_line() +
    scale_x_continuous("Day", breaks = seq(0, 25, 5)) +
    scale_y_continuous("Pr(recovery)") +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    ) 
  
  l_res_1$rmsts$trt <- factor(l_res_1$rmsts$trt, levels = c("soc", "defer", "discont"))
  p2 <- ggplot(l_res_1$rmsts, aes(x = trt, y = rmst_mu)) +
    geom_point() +
    geom_linerange(aes(ymin = rmst_lo, ymax = rmst_hi)) +
    scale_x_discrete("") +
    scale_y_continuous("RMST") +
    theme_minimal() 
   
  p1 + p2
  

  d_res_2_defer <- compare_rmst(
    trt_a = "soc", trt_b = "defer",
    
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    # episode window for recovery metrics
    t_window   = l_spec$rmst_eh_horizon,
    delta_ni = 0.75,
    # sample data for covariate profile
    d_cohort
  )
  d_res_2_discont <- compare_rmst(
    trt_a = "soc", trt_b = "discont",
    
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    # episode window for recovery metrics
    t_window   = l_spec$rmst_eh_horizon,
    delta_ni = 0.75,
    # sample data for covariate profile
    d_cohort
  )
  d_res_2 <- rbind(d_res_2_defer, d_res_2_discont)
  kableExtra::kbl(d_res_2, format = "simple")
  d_res_2[, contrast := paste0(trt_b, " - ", trt_a)]
  p1 <- ggplot(d_res_2, aes(x = contrast, y = delta_mu)) +
    geom_point() +
    geom_linerange(aes(ymin = delta_lo, ymax = delta_hi)) +
    scale_x_discrete("") +
    scale_y_continuous("Difference in RMST") +
    theme_minimal() 
  p1
  #
  
  ix_sample <- sample(1:length(po_he_shape), size = 500, replace = F)
  # ix_sample <- 1:length(po_he_shape)
  tictoc::tic()
  l_res_3 <- compare_trajectory(
    trt_a = "soc",
    trt_b = "defer",
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_defer, po_eh_b_discont, po_eh_u_a,
    ix_sample,
    
    followup  = 365,
    max_trans = 50L,
    d_cohort
  )
  tictoc::toc()
  l_res_3
  
  
  
  
  
  
}



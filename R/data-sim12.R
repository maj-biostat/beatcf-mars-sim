suppressWarnings(library(data.table))
suppressWarnings(library(mvtnorm))
# for weibullph
suppressWarnings(library(flexsurv))
suppressWarnings(library(ggplot2))
suppressWarnings(library(INLAjoint))
suppressWarnings(library(pbapply))
suppressWarnings(library(lobstr))

suppressWarnings(library(patchwork))

# todo 
# review betancourt warping paper etc



# Baseline ppfev with option for subject level heterogeneity
# ppfev_0 <- function(age, A = 123, k = 0.4, p = 0.3, sd_ppfev = 3) {
#   A * exp(-k * age^p) + rnorm(length(age), 0, sd_ppfev)
# }

# Baseline ppfev with option for subject level heterogeneity
ppfev_0 <- function(age, age_min = 10, ppfev0_max = 100, k = 0.1, p = 0.5, sd_ppfev = 0) {
  ppfev0_max * exp(-k * (age-age_min)^p) + rnorm(length(age), 0, sd_ppfev)
}




# age_mean_log = log(35)
# age_sd_log = 0.4
# age_min = 10
# age_max = 60
# # frailty parameters to link recurrences
# sd_he = 0.3
# sd_eh = 0.3 
# rho_frailty = -0.4
# shape_he = 1.1 
# shape_eh = 0.9
# ppfev_ref = 77.5
# # linear predictor exacerbation
# mu_exacerb = -4.5 
# b_ppfev_exacerb = -0.02
# # linear predictor recovery
# mu_recov = -0.5 
# b_ppfev_recov = 0.01
# # trt is soc, delay, defer
# b_trt = c(0, -0.3, -0.2)

get_sim12_pt <- function(
    l_spec
  ){
  
  # baseline 
  # probably should include sex in this...?
  age <- rlnorm(1, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd)
  age <- pmin(pmax(age, l_spec$age_min), l_spec$age_max)
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


get_sim12_stan_data <- function(d_all){
  
  
  
  # The analysis proceeds on only those that have exacerbations
  d_tmp <- copy(d_all[state == "E"])
  d_tmp[, ppfev_std := ppfev_baseline]
  d_tmp[, delay := as.numeric(trt == "delay")]
  d_tmp[, defer := as.numeric(trt == "defer")]
  
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
  
  
  X_obs <- model.matrix(~-1 + ppfev_std + delay + defer, data = d_obs)
  X_cens <- model.matrix(~-1 + ppfev_std + delay + defer, data = d_cens)
  
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
    pri_s_u = 1
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
    pri_s_u = 1
  )
  
  list(
    ld_eh = ld_eh,
    ld_he = ld_he
  )
  
}

sim_study <- function(){
  
  n_sim <- 100
  smry_list <- list()
  
  l_spec <- get_demo_spec()
  
  i <- 1
  while(i <= n_sim){
    
    d_cohort <- get_sim12_cohort(l_spec)
  
    # for INLAjoint, I believe data.table causes some challenges
    d <- data.frame(
      id = d_cohort$id,
      gap_time = d_cohort$dur,
      evt_he = as.numeric(d_cohort$state == "H"),
      evt_eh = as.numeric(d_cohort$state == "E"),
      ppfev_std = d_cohort$ppfev_baseline,
      trt = factor(d_cohort$trt, level = c("soc", "defer", "delay", "none"))
      
    )
    
    f1 <- joint(
      formSurv =
        list(
          inla.surv(time = gap_time, event = evt_he) ~ ppfev_std + (1 | id),
          inla.surv(time = gap_time, event = evt_eh) ~ ppfev_std + trt + (1 | id)
        ),
      basRisk= "weibullsurv", id = "id", dataSurv = d, control= list(config=TRUE)
    )
    
    d_smry <- rbind(
      cbind(
        model = "S1",
        parameter = rownames(summary(f1)$SurvEff[[1]]),
        data.table(summary(f1)$SurvEff[[1]], hr = T)
      ),
      cbind(
        model = "S1",
        parameter = rownames(summary(f1)$ReffListS[[1]]),
        data.table(summary(f1)$ReffListS[[1]])
      ),
      cbind(
        model = "S2",
        parameter = rownames(summary(f1)$SurvEff[[2]]),
        data.table(summary(f1)$SurvEff[[2]], hr = T)
      ),
      cbind(
        model = "S2",
        parameter = rownames(summary(f1)$ReffListS[[2]]),
        data.table(summary(f1)$ReffListS[[2]])
      ),
      fill = T
    )
    
    smry_list[[i]] <- d_smry
    
    i <- i + 1
    
  }
  
  d_smry <- rbindlist(smry_list, idcol = "id_sim")
  
  # just extract means and see if they have any resemblence to the dgp params
  d_smry[, .(
    mu = mean(mean), 
    q_025 = mean(`0.025quant`),
    q_975 = mean(`0.975quant`),
    mu_min = min(mean),
    mu_max = max(mean)
  ), keyby = .(model, parameter)]
 
}


example_inla <- function(){
  
  
  library(p3state.msm)
  data(heart2)
  
  head(heart2, 10)
  
  # pt on waiting list and can be transplated and then die or just one or none
  # (continue in waiting list)
  
  # times1: time of transplant/censoring time.
  # delta: transplant indicator (1:yes; 0:no).
  # times2: time to death since the transplant/censoring time.
  # time: times1+times2.
  # status: censoring indicator (1:dead; 0:alive).
  # age: age-48 years.
  # year: year of acceptance (in years after 1Nov1967).
  # surgery: prior bypass surgery (1:yes; 0:no).
  
  # first col indicates transplant event
  # second col indicates waiting list to death without transplant
  # third col indicates death following transplant
  # row with all zeros indicates censoring while on waiting list
  event <- matrix(c(
    heart2$delta, 
    heart2$status * (1 - heart2$delta),
    heart2$status * heart2$delta 
    ), 
    ncol= 3)
  head(event)
  
  heart2$times3 <- ifelse(heart2$times2== 0, heart2$times1, heart2$times2)
  
  
  m5.ms <- joint(formSurv= list(
    inla.surv(times1, event[,1]) ~ age + year + surgery,
    inla.surv(times3, event[,2]) ~ age + year + surgery,
    inla.surv(times2, event[,3])~ age + year + surgery),
    basRisk= rep("weibullsurv", 3), dataSurv= heart2)
  
  
  summary(m5.ms)
  
  # surv.obj.ms3 <- inla.surv(time= heart2$time, truncation= heart2$times1,
  #                           event = event[,3])
  
  
  
  t <- seq(0.1, 1000, by= 1)
  riskW <- function(t, lambda, alpha) {
    # instantaneous risk under weibullPH parameterisation
    # with shape a and scale lambda
    alpha*lambda*t^(alpha-1)
  }
  risk1 <- riskW(t, 
                 lambda = exp(m5.ms$summary.fixed["Intercept_S1", "mean"]),
                 alpha = m5.ms$summary.hyperpar$mean[1])
  risk2 <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S2", "mean"]),
                 m5.ms$summary.hyperpar$mean[2])
  
  risk3 <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S3", "mean"]),
                 m5.ms$summary.hyperpar$mean[3])
  
  risk1s <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S1", "mean"] +
                           m5.ms$summary.fixed["surgery_S1", "mean"]),
                  m5.ms$summary.hyperpar$mean[1])
  risk2s <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S2", "mean"] +
                           m5.ms$summary.fixed["surgery_S2", "mean"]),
                  m5.ms$summary.hyperpar$mean[2])
  risk3s <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S3", "mean"] +
                           m5.ms$summary.fixed["surgery_S3", "mean"]),
                  m5.ms$summary.hyperpar$mean[3])
  
  # probability of staying in state 1 given currently in state 1
  # transition rate matrix is of the form
  #     | -(\rho_{12} + \rho_{13})   \rho_{12}                  \rho_{13}                |
  # Q = | \rho_{21}                -(\rho_{21}  + \rho_{23})    \rho_{23}                |
  #     | \rho_{31}                  \rho_{32}                -(\rho_{31}  + \rho_{32})  |
  
  # state transition from exp(Q \Delta) where \Delta is the time interval
  # here the increment is 1 unit
  
  # remaining on waiting list
  p11 <- exp(-cumsum(risk1)-cumsum(risk2))
  p11s <- exp(-cumsum(risk1s)-cumsum(risk2s))
  
  # remaining in the transplant state
  p22 <- exp(-cumsum(risk3))
  p22s <- exp(-cumsum(risk3s))
  
  # transition from waiting list to transplant
  p12 <- cumsum(p11*risk1*p22)
  plot(t, p12, type ="l")
  p12s <- cumsum(p11s*risk1s*p22s)
  
  # transition from waiting list to death
  p13 <- 1-p11-p12
  p13s <- 1-p11s-p12s
  
  # transition from transplant to death
  p23 <- 1-p22
  p23s <- 1-p22s
  #
  
  
  
}


example_hesim <- function(){
  
  library("hesim")
  library("data.table")
  
  # Treatment strategies
  strategies <- data.table(
    strategy_id = c(1, 2),
    strategy_name  = c("SOC", "New"))
  
  # Patients
  n_patients <- 1000
  patients <- data.table(patient_id = 1:n_patients,
                         age = rnorm(n_patients, mean = 45, sd = 7),
                         female = rbinom(n_patients, size = 1, prob = .51))
  patients[, grp_id := ifelse(female == 1, 1, 2)]
  patients[, grp_name := ifelse(female == 1, "Female", "Male")]
  
  # (Non-death) health states
  states <- data.table(state_id = c(1, 2),
                       state_name = c("Stage 1", "Stage 2")) 
  
  # Transitions
  tmat <- rbind(c(NA, 1, 2),
                c(3, NA, 4),
                c(NA, NA, NA))
  colnames(tmat) <- rownames(tmat) <- c("Stage 1", "Stage 2", "Death")
  transitions <- create_trans_dt(tmat)
  transitions[, trans := factor(transition_id)]
  
  # Combining
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients, 
                          states = states,
                          transitions = transitions)
  print(hesim_dat)
  
  labs <- get_labels(hesim_dat)
  print(labs)
  
  
  library("flexsurv")
  mstate_data <- data.table(mstate3_exdata$transitions)
  mstate_data[, trans := factor(trans)]
  mstate_data[1:20,][order(patient_id, Tstart)]
  fit_wei <- flexsurv::flexsurvreg(Surv(years, status) ~ trans + 
                                     factor(strategy_id):trans +
                                     age:trans + 
                                     female: trans +
                                     shape(trans), 
                                   data = mstate_data, 
                                   dist = "weibull")
  #
  # Utility
  utility_tbl <- stateval_tbl(
    data.table(state_id = states$state_id,
               mean = mstate3_exdata$utility$mean,
               se = mstate3_exdata$utility$se),
    dist = "beta"
  )
  
  # Costs
  drugcost_tbl <- stateval_tbl(
    data.table(strategy_id = strategies$strategy_id,
               est = mstate3_exdata$costs$drugs$costs),
    dist = "fixed"
  )
  
  medcost_tbl <- stateval_tbl(
    data.table(state_id = states$state_id,
               mean = mstate3_exdata$costs$medical$mean,
               se = mstate3_exdata$costs$medical$se),
    dist = "gamma"
  )
  
  n_samples <- 1000
  transmod_data <- expand(hesim_dat, 
                          by = c("strategies", "patients", "transitions"))
  head(transmod_data)
  
  transmod <- create_IndivCtstmTrans(fit_wei, transmod_data,
                                     trans_mat = tmat, n = n_samples,
                                     uncertainty = "normal")
  class(transmod)
  
  
  utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)
  
  # Costs
  drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat)
  medcostmod <- create_StateVals(medcost_tbl, n = n_samples, hesim_data = hesim_dat)
  costmods <- list(drugs = drugcostmod,
                   medical = medcostmod)
  
  ictstm <- IndivCtstm$new(trans_model = transmod,
                           utility_model = utilitymod,
                           cost_models = costmods)
  
  # creates 1000 samples of 1000 pt trajectories
  ictstm$sim_disease()
  ictstm$disprog_[1:20]
  unique(ictstm$disprog_$sample)
  unique(ictstm$disprog_$patient_id)
  
  d_fig <- ictstm$disprog_[sample == 1 & patient_id %in% 1:10]
  
  library(ggplot2)  
  ggplot(d_fig, aes(x = time_stop, y = to)) +
    geom_point() + 
    geom_line() +
    facet_wrap(~patient_id) 
}






simulate_occupancy <- function(
    # posterior draws
  po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
  po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
  po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
  # simulation settings
  S         = 100,
  N_pop     = 1000,
  followup  = 365,
  eval_times = seq(0, 365, by = 1),
  # covariate profile
  ppfev_std = 0,
  trt       = "soc",   # "soc", "delay", "defer"
  # max transitions per patient (safety cap)
  max_trans = 200L
) {
  
  n_eval  <- length(eval_times)
  n_draws <- length(po_he_shape)
  
  # sample S draw indices without replacement (or with, if S > n_draws)
  draw_idx <- sample(n_draws, S, replace = S > n_draws)
  
  # treatment coefficient selector (applied to EH transition)
  b_trt_val <- switch(trt,
                      soc   = 0,
                      delay = 1,  # will be multiplied by b_delay below
                      defer = 1   # will be multiplied by b_defer below
  )
  
  # result matrix: rows = eval_times, cols = posterior draws
  # store proportion in H at each eval time for each draw
  occ_H <- matrix(NA_real_, nrow = n_eval, ncol = S)
  
  s <- 1
  for (s in seq_len(S)) {
    
    d <- draw_idx[s]
    
    he_shape <- po_he_shape[d]
    he_b0    <- po_he_b_0[d]
    he_bfev  <- po_he_b_ppfev[d]
    he_ua    <- po_he_u_a[d]
    
    eh_shape <- po_eh_shape[d]
    eh_b0    <- po_eh_b_0[d]
    eh_bfev  <- po_eh_b_ppfev[d]
    eh_btrt  <- switch(trt,
                       soc   = 0,
                       delay = po_eh_b_delay[d],
                       defer = po_eh_b_defer[d]
    )
    eh_ua <- po_eh_u_a[d]
    
    u_he <- rgamma(N_pop, shape = he_ua, rate = he_ua)
    u_eh <- rgamma(N_pop, shape = eh_ua, rate = eh_ua)
    
    # at given ppfev_std (scalar here, could be a vector for a population)
    lam_he <- u_he * exp(he_b0 + he_bfev * ppfev_std)
    lam_eh <- u_eh * exp(eh_b0 + eh_bfev * ppfev_std + eh_btrt)
    
    # Simulate enough transitions to cover followup for all patients.
    # Strategy: simulate max_trans sojourn times per patient in bulk,
    # compute cumulative transition times, then look up state at eval_times.
    
    # HE transitions: sojourn in H (transitions 1, 3, 5, ...)  -> N_pop x max_trans/2
    # EH transitions: sojourn in E (transitions 2, 4, 6, ...)  -> N_pop x max_trans/2
    n_he <- ceiling(max_trans / 2)
    n_eh <- floor(max_trans / 2)
    
    # Simulate all sojourn times at once: N_pop x n_he matrices
    # rweibullPH(n, shape, scale): scale = lambda in PH form
    # Equivalent: (-log(U) / lambda)^(1/shape)
    
    # HE sojourns — one row per patient, one col per transition
    U_he <- matrix(runif(N_pop * n_he), nrow = N_pop)
    # broadcasting lam_he as column
    soj_he <- (-log(U_he) / lam_he)^(1 / he_shape)   
    
    U_eh <- matrix(runif(N_pop * n_eh), nrow = N_pop)
    soj_eh <- (-log(U_eh) / lam_eh)^(1 / eh_shape)
    
    # Interleave HE and EH columns: H1, E1, H2, E2, ...
    # Transition times are cumulative sums of interleaved sojourns
    # Interleaved sojourn matrix: N_pop x max_trans
    soj_interleaved <- matrix(NA_real_, nrow = N_pop, ncol = n_he + n_eh)
    he_cols <- seq(1, by = 2, length.out = n_he)
    eh_cols <- seq(2, by = 2, length.out = n_eh)
    soj_interleaved[, he_cols] <- soj_he
    soj_interleaved[, eh_cols] <- soj_eh
    
    # Cumulative transition times: N_pop x max_trans
    cum_times <- t(apply(soj_interleaved, 1, cumsum))
    
    # State at each transition boundary alternates H->E->H->E...
    # Col k of cum_times is when the patient *leaves* state:
    #   col 1 (he): leaves H, enters E
    #   col 2 (eh): leaves E, enters H
    #   etc.
    # State just AFTER transition k:
    #   odd k -> entered E
    #   even k -> entered H
    # So state at time t: count how many transitions have occurred by t (= n_crossed),
    # if n_crossed is even -> in H, if odd -> in E
    
    # For each patient and eval time, count transitions crossed
    # eval_times: length n_eval
    # cum_times:  N_pop x max_trans
    # Result: N_pop x n_eval matrix of states
    
    in_H <- matrix(NA, nrow = N_pop, ncol = n_eval)
    
    for (e in seq_len(n_eval)) {
      tt <- eval_times[e]
      # number of transitions crossed by time tt, per patient
      n_crossed <- rowSums(cum_times <= tt)
      # even crossings -> in H (started in H)
      in_H[, e] <- (n_crossed %% 2L == 0L)
    }
    
    occ_H[, s] <- colMeans(in_H)
    
  }
  
  # Summarise across posterior draws
  data.table(
    time     = eval_times,
    mean_H   = rowMeans(occ_H),
    lo_H     = apply(occ_H, 1, quantile, 0.025),
    hi_H     = apply(occ_H, 1, quantile, 0.975),
    mean_E   = 1 - rowMeans(occ_H),
    lo_E     = 1 - apply(occ_H, 1, quantile, 0.975),
    hi_E     = 1 - apply(occ_H, 1, quantile, 0.025)
  )
}

# for each posterior draw, obtain the durations times for exacerbation for a
# population of patients characterised by the covariate specification and
# nominated treatment assignment
# compute the proportion recovered over the evaluation interval to the 
# max followup period stipulated in t_window
# the proportion recovered is our proxy for the probability of recovered
# with uncertainty propagated from the posterior draws.
# compute the rmst (to t_window) as the average of the exacerbation duration
simulate_episode <- function(
    # posterior draws (EH model only needed for primary metric)
  po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
  po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
  # simulation settings
  S          = 200,
  N_pop      = 2000,
  # episode window for recovery metrics 
  # (for bound on restricted mean time to recovery)
  t_window   = 30,
  # evaluate at each time increment
  eval_times = seq(0, 30, by = 1),
  trt        = "soc",
  # sample data
  d_cohort,
  m_ix
  
) {
  
  n_eval  <- length(eval_times)
  n_draws <- length(po_eh_shape)
  draw_idx <- sample(n_draws, S, replace = S > n_draws)
  
  # storage: recovery probability at each eval time, per draw
  # and RMST per draw
  p_recovered <- matrix(NA_real_, nrow = n_eval, ncol = S)
  rmst        <- numeric(S)
  
  s <- 1
  for (s in seq_len(S)) {
    
    d <- draw_idx[s]
    
    eh_shape <- po_eh_shape[d]
    eh_ua    <- po_eh_u_a[d]
    eh_btrt  <- switch(trt,
                       soc   = 0,
                       delay = po_eh_b_delay[d],
                       defer = po_eh_b_defer[d]
    )
    
    # scale parameter: mu = exp(b0 + b_ppfev * ppfev_std + b_trt)
    log_mu <- po_eh_b_0[d] + po_eh_b_ppfev[d] * d_cohort$ppfev_baseline[m_ix[,s]] + eh_btrt
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


# simulates transitions H -> E and E -> H for a population of patients in order 
# to estimate the total number of exacerbations and the mean total duration of 
# exacerbations over the year. uncertainty is propogated from the posterior 
# draws
simulate_trajectory <- function(
    # HE posterior draws
  po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
  # EH posterior draws
  po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
  po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
  # settings
  S         = 200,
  N_pop     = 2000,
  followup  = 365,
  trt       = "soc",
  max_trans = 300L,
  draw_idx_override = NULL,
  d_cohort,
  m_ix
) {
  
  if(is.null(draw_idx_override)){
    n_draws  <- length(po_he_shape)
    draw_idx <- sample(n_draws, S, replace = S > n_draws)
  } else {
    n_draws  <- length(draw_idx_override)
    draw_idx <- draw_idx_override
    # reset S based on fixed set of draws
    S <- length(draw_idx)
  }
  
  
  # per-draw summaries
  mean_time_E  <- numeric(S)
  mean_n_exac  <- numeric(S)
  
  s <- 1
  for (s in seq_len(S)) {
    
    d <- draw_idx[s]
    
    # parameters
    he_shape <- po_he_shape[d]
    he_ua    <- po_he_u_a[d]
    he_mu    <- exp(po_he_b_0[d] + po_he_b_ppfev[d] * d_cohort$ppfev_baseline[m_ix[, s]])
    
    eh_btrt  <- switch(trt,
                       soc   = 0,
                       delay = po_eh_b_delay[d],
                       defer = po_eh_b_defer[d]
    )
    eh_shape <- po_eh_shape[d]
    eh_ua    <- po_eh_u_a[d]
    eh_mu    <- exp(po_eh_b_0[d] + po_eh_b_ppfev[d] * d_cohort$ppfev_baseline[m_ix[, s]] + eh_btrt)
    
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

# contrast for mean time in E state and mean number of exacerbation


# Call episodes once per arm, then contrast
compare_treatments <- function(
    # arms to compare
    arms = c("soc", "delay", "defer"), 
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    # simulation settings
    S          = 200,
    N_pop      = 2000,
    # episode window for recovery metrics
    t_window   = 30,
    eval_times = seq(0, 30, by = 1),
    # sample data
    d_cohort
    ) {
  
  # cohort covariate indexes to use in prediction
  m_ix <- sapply(1:S, function(ii){
    sample(1:nrow(d_cohort), size = N_pop, replace = T)  
  })
  
  trt <- arms[1]
  results <- lapply(arms, function(trt) {
    
    res <- simulate_episode(
      po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
      po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
      # simulation settings
      S , N_pop ,
      t_window, eval_times,
      # simulate each arm independently
      trt = trt,
      # sample data
      d_cohort,
      # which rows from the cohort data should be used in prediction
      m_ix
      )
    res$recovery_curve[, trt := trt]
    list(curve = res$recovery_curve, rmst = cbind(trt = trt, res$rmst))
  })
  
  curves <- rbindlist(lapply(results, `[[`, "curve"))
  rmsts  <- rbindlist(lapply(results, `[[`, "rmst"))
  
  # pairwise RMST contrasts (delay vs soc, defer vs soc)
  # re-run with paired draws to get proper posterior contrast
  list(curves = curves, rmsts = rmsts)
}

# compares the rmst to the specified t_window limit across treatment groups
# under assumed covariate pattern
# estimates obtained by repeatedly simulating sourjoun times for the 
# exacerbation state for a population of patients under each draw of the 
# posterior with the same draws being used for each treatment group.
contrast_rmst <- function(
    trt_a = "soc", trt_b = "delay",
    
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    # simulation settings
    S          = 200,
    N_pop      = 2000,
    # episode window for recovery metrics
    t_window   = 30,
    delta_ni = 0,
    # sample data for covariate profile
    d_cohort
) {
  
  n_draws  <- length(po_eh_shape)
  draw_idx <- sample(n_draws, S, replace = S > n_draws)
  
  # cohort covariate indexes to use in prediction
  m_ix <- sapply(1:S, function(ii){
    sample(1:nrow(d_cohort), size = N_pop, replace = T)  
  })
  
  
  rmst_a <- rmst_b <- numeric(S)
  
  for (s in seq_len(S)) {
    d <- draw_idx[s]
    
    for (trt in c(trt_a, trt_b)) {
      eh_btrt <- switch(trt, soc = 0,
                        delay = po_eh_b_delay[d],
                        defer = po_eh_b_defer[d])
      mu  <- exp(po_eh_b_0[d] + po_eh_b_ppfev[d] * d_cohort$ppfev_baseline[m_ix[, s]] + eh_btrt)
      u_i <- rgamma(N_pop, po_eh_u_a[d], po_eh_u_a[d])
      soj <- (-log(runif(N_pop)) / (u_i * mu))^(1 / po_eh_shape[d])
      val <- mean(pmin(soj, t_window))
      if (trt == trt_a) rmst_a[s] <- val else rmst_b[s] <- val
    }
  }
  
  contrast <- rmst_b - rmst_a
  data.table(
    trt_a    = trt_a,
    trt_b    = trt_b,
    delta_mu     = mean(contrast),
    delta_lo       = quantile(contrast, 0.025),
    delta_hi       = quantile(contrast, 0.975),
    # does trt_b result in a rmst that is above the rmst for trt_a by less
    # than the ni margin
    pr_lt_ni = mean(contrast < delta_ni)   
  )
}


contrast_trajectory <- function(
    trt_a = "soc",
    trt_b = "delay",
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    # settings
    S         = 200,
    N_pop     = 2000,
    followup  = 365,
    max_trans = 300L,
    d_cohort
) {
  # pass same draw_idx to both by fixing the RNG seed approach:
  # easier to just run both inside one loop with shared draw_idx
  
  n_draws  <- length(po_he_shape)
  draw_idx <- sample(n_draws, S, replace = S > n_draws)
  
  
  
  # cohort covariate indexes to use in prediction
  m_ix <- sapply(1:S, function(ii){
    sample(1:nrow(d_cohort), size = N_pop, replace = T)  
  })
  
  # po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
  # # EH posterior draws
  # po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
  # po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
  # settings
  # S         = 200
  # N_pop     = 2000
  # followup  = 365
  # ppfev_std = 0
  # trt       = "soc"
  # max_trans = 300L
  # draw_idx_override = NULL
  
  res_a <- simulate_trajectory(
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    S, N_pop, followup, 
    trt = trt_a,
    max_trans,
    draw_idx_override = draw_idx,
    d_cohort,
    m_ix = m_ix
  )
  
  res_b <- simulate_trajectory(
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    S ,N_pop, followup, 
    trt = trt_b,
    max_trans,
    draw_idx_override = draw_idx,
    d_cohort,
    m_ix = m_ix
  )
  
  contrast_time  <- res_b$draws$mean_time_E - res_a$draws$mean_time_E
  contrast_nexac <- res_b$draws$mean_n_exac - res_a$draws$mean_n_exac
  
  list(
    time_E = data.table(
      trt_a       = trt_a,  
      trt_b = trt_b,
      mean_E_trt_a = res_a$time_E$mean,
      mean_E_trt_b = res_b$time_E$mean,
      delta_mu        = mean(contrast_time),
      delta_lo          = quantile(contrast_time, 0.025),
      delta_hi          = quantile(contrast_time, 0.975),
      pr_positive = mean(contrast_time > 0)
    ),
    n_exac = data.table(
      trt_a       = trt_a,  trt_b = trt_b,
      n_exac_trt_a = res_a$n_exac$mean,
      n_exac_trt_b = res_b$n_exac$mean,
      delta_mu        = mean(contrast_nexac),
      delta_lo          = quantile(contrast_nexac, 0.025),
      delta_hi          = quantile(contrast_nexac, 0.975),
      pr_positive = mean(contrast_nexac > 0)
    )
  )
}



weibullPH_stats <- function(shape_ph, scale_ph, tau = 25){
  # aft scale
  aft_lambda <- scale_ph^(-1/shape_ph)
  # aft shape
  aft_k <- shape_ph
  # median
  w_med <- aft_lambda * (log(2))^{1/aft_k}
  # mean
  w_mu <- aft_lambda * gamma(1 + (1/aft_k))
  # variance
  w_sd <- sqrt((aft_lambda^2)*(gamma(1 + 2/aft_k) - (gamma(1 + 1/aft_k))^2))
  
  SweibullPH <- function(x, shape_ph, scale_ph){
    1-flexsurv::pweibullPH(x, shape_ph, scale_ph)
  }
  
  # plot(0:365, SweibullPH(0:365, shape_ph, scale_ph))
  
  rmst <- integrate(SweibullPH, 0, tau, shape_ph = shape_ph, scale_ph = scale_ph)
  
  
  c(median = w_med, mean = w_mu, sd = w_sd, rmst = rmst$value)
}

example_stan_2 <- function(){
  
  
  m3 <- cmdstanr::cmdstan_model("stan/sim12-v03.stan")
  
  l_spec <- get_demo_spec()
  
  l_spec$shape_he <- 2.85 # originally 1.1
  l_spec$mu_exacerb <- -16.1 # originally  -4.5
  weibullPH_stats(
    shape_ph = l_spec$shape_he, scale_ph = exp(l_spec$mu_exacerb), tau = 365)

  y_he <- flexsurv::rweibullPH(
    1e4, l_spec$shape_he, scale = exp(l_spec$mu_exacerb)
  )
  hist(y_he, xlim = c(0, max(y_he)*1.1))
  c(median(y_he), mean(y_he), sd(y_he))
  
  l_spec$shape_eh <- 2.75  # originally 0.9 
  l_spec$mu_recov <- -6.5  # originally -0.5
  (res_1 <- weibullPH_stats(
    shape_ph = l_spec$shape_eh, scale_ph = exp(l_spec$mu_recov), tau = 25))
  
  log_hr <- seq(0.12,0.15, len = 100)
  delta <- numeric(100)
  for(ii in seq_along(log_hr)){
    (res_2 <- weibullPH_stats(
      shape_ph = l_spec$shape_eh, scale_ph = exp(l_spec$mu_recov - log_hr[ii]), tau = 25))
    delta[ii] <- res_2["rmst"] - res_1["rmst"]
  }
  names(delta) <- paste0(log_hr)
  delta
  
  (res_2 <- weibullPH_stats(
    shape_ph = l_spec$shape_eh, scale_ph = exp(l_spec$mu_recov - 0.14168), tau = 25))
  res_2["rmst"] - res_1["rmst"]
  
  y_eh <- flexsurv::rweibullPH(
    1e4, l_spec$shape_eh, scale = exp(l_spec$mu_recov)
  )
  hist(y_eh)
  #
  
  # some of the cohort will have no exacerbation and so the event indicator 
  # for the healthy state record should have an event indicator set to zero 
  # due to administrative censoring
  
  # similarly, the last record for each pt will almost certainly be censored
  # again, this is administrative censoring at 365 days
  
  
  d_cohort <- get_sim12_cohort(l_spec)
  
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
  d_tmp[, delay := as.numeric(trt == "delay")]
  d_tmp[, defer := as.numeric(trt == "defer")]
  
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
  
  
  X_obs <- model.matrix(~-1 + ppfev_std + delay + defer, data = d_obs)
  X_cens <- model.matrix(~-1 + ppfev_std + delay + defer, data = d_cens)
  
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
    pri_s_u = 1
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
    pri_s_u = 1
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
  f_1_1$summary(variables = c("shape", "b_0", "b", "u_a"))
  
  f_1_2 <- m3$sample(
    ld_2, iter_warmup = 500, iter_sampling = 500,
    parallel_chains = 1, chains = 1, refresh = 10, show_exceptions = T,
    max_treedepth = 11,
    output_dir = output_dir_mcmc,
    output_basename = foutname
  )
  f_1_2$summary(variables = c("shape", "b_0", "b", "u_a"))
  # c(l_spec$shape_he, l_spec$shape_eh)
  # c(l_spec$mu_exacerb, l_spec$mu_recov)
  # c(l_spec$b_trt)
  # l_spec$b_ppfev_recov
  
  # laplace approx
  f_2_1_optim <- m3$optimize(data = ld_1, jacobian = TRUE)
  f_2_1 <- m3$laplace(data = ld_1, mode = f_2_1_optim, draws = 2000)
  f_2_1$summary(variables = c("shape", "b_0", "b", "u_a"))
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
  po_eh_b_delay <- d_post_eh$`b[2]`
  po_eh_b_defer <- d_post_eh$`b[3`
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
  
  
  l_res_1 <- compare_treatments(
    arms = c("soc", "delay", "defer"),
    po_eh_shape, 
    po_eh_b_0, po_eh_b_ppfev, po_eh_b_delay, po_eh_b_defer, 
    po_eh_u_a,
    # simulation settings
    S          = 200,
    N_pop      = 2000,
    # episode window for recovery metrics
    t_window   = l_spec$rmst_eh_horizon,
    eval_times = seq(0, l_spec$rmst_eh_horizon, by = 1),
    # sample data for covariate profile
    d_cohort
  )
  l_res_1
  
  library(ggplot2)
  p1 <- ggplot(l_res_1$curves, aes(x = time, y = mean, group = trt, col = trt)) +
    geom_line() +
    scale_x_continuous("Day", breaks = seq(0, 25, 5)) +
    scale_y_continuous("Pr(recovery)") +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    ) 
  
  l_res_1$rmsts$trt <- factor(l_res_1$rmsts$trt, levels = c("soc", "defer", "delay"))
  p2 <- ggplot(l_res_1$rmsts, aes(x = trt, y = rmst_mu)) +
    geom_point() +
    geom_linerange(aes(ymin = rmst_lo, ymax = rmst_hi)) +
    scale_x_discrete("") +
    scale_y_continuous("RMST") +
    theme_minimal() 
   
  p1 + p2
  
  
  

  d_res_2_delay <- contrast_rmst(
    trt_a = "soc", trt_b = "delay",
    
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    # simulation settings
    S          = 200,
    N_pop      = 2000,
    # episode window for recovery metrics
    t_window   = l_spec$rmst_eh_horizon,
    delta_ni = 0.75,
    # sample data for covariate profile
    d_cohort
  )
  d_res_2_defer <- contrast_rmst(
    trt_a = "soc", trt_b = "defer",
    
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    # simulation settings
    S          = 200,
    N_pop      = 2000,
    # episode window for recovery metrics
    t_window   = l_spec$rmst_eh_horizon,
    delta_ni = 0.75,
    # sample data for covariate profile
    d_cohort
  )
  d_res_2 <- rbind(d_res_2_delay, d_res_2_defer)
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
  
  
  #####
  l_res_3 <- contrast_trajectory(
    trt_a = "soc",
    trt_b = "delay",
    po_he_shape, po_he_b_0, po_he_b_ppfev, po_he_u_a,
    # EH posterior draws
    po_eh_shape, po_eh_b_0, po_eh_b_ppfev,
    po_eh_b_delay, po_eh_b_defer, po_eh_u_a,
    # simulation settings
    S          = 200,
    N_pop      = 2000,
    followup  = 365,
    max_trans = 300L,
    d_cohort
  )
  l_res_3
  
  
  
  
  
  
}

example_stan_1 <- function(){
  
  
  m1 <- cmdstanr::cmdstan_model("stan/sim12-v02.stan")
  
  l_spec <- get_demo_spec()
  l_spec$N_pt <- 1000
  l_spec$is <- 1
  l_spec$ie <- 1000
  
  
  d_cohort <- get_sim12_cohort(l_spec)
  
  # 
  # int<lower=0> N_obs_he;
  # int N_id_he;
  # vector[N_obs_he] y_obs_he;
  # array[N_obs_he] int id_he;
  # // ppfev (mean centred)
  # vector[N_obs_he] ppfev_he;
  # 
  # int<lower=0> N_obs_eh;
  # int N_id_eh;
  # vector[N_obs_eh] y_obs_eh;
  # array[N_obs_eh] int id_eh;
  # // ppfev (mean centred)
  # vector[N_obs_eh] ppfev_eh;
  # array[N_obs_eh] int trt_ix_eh;
  
  # d_tmp <- d_cohort[state == "E"]
  # ld <- list(
  #   N_obs_he = d_cohort[state == "H", .N],
  #   N_id_he = d_cohort[state == "H", length(unique(id))],
  #   y_obs_he = d_cohort[state == "H", dur],
  #   id_he = d_cohort[state == "H", id],
  #   ppfev_he = d_cohort[state == "H", ppfev_baseline],
  #   
  #   N_obs_eh = d_cohort[state == "E", .N],
  #   N_id_eh = d_cohort[state == "E", length(unique(id))],
  #   y_obs_eh = d_cohort[state == "E", dur],
  #   id_eh = d_cohort[state == "E", id],
  #   ppfev_eh = d_cohort[state == "E", ppfev_baseline],
  #   trt_ix_eh = d_cohort[state == "E", factor(trt, levels = c("soc", "delay", "defer"))]
  # )
  
  
  # int<lower=0> N_obs;
  # int N_id;
  # vector[N_obs] y_obs;
  # array[N_obs] int id;
  # // ppfev (mean centred)
  # vector[N_obs] ppfev;
  # 
  d_tmp <- copy(d_cohort[state == "E"])
  d_tmp[, ppfev_std := ppfev_baseline]
  d_tmp[, delay := as.numeric(trt == "delay")]
  d_tmp[, defer := as.numeric(trt == "defer")]
  d_tmp[, id_idx := as.integer(factor(id))]
  X <- model.matrix(
    ~-1 + ppfev_std + delay + defer,
    data = d_tmp)
  
  ld_1 <- list(
    N_obs = nrow(d_tmp),
    N_id = length(unique(d_tmp$id)),
    y_obs = d_tmp$dur,
    P = ncol(X),
    X = X,
    # ensure ids are contiguous (some pt might not have pex)
    id = d_tmp$id_idx,
    pri_s_u = 1
  )
  
  d_tmp <- copy(d_cohort[state == "H"])
  d_tmp[, ppfev_std := ppfev_baseline]
  d_tmp[, id_idx := as.integer(factor(id))]
  X <- model.matrix(
    ~-1 + ppfev_std,
    data = d_tmp)
  
  ld_2 <- list(
    N_obs = nrow(d_tmp),
    N_id = length(unique(d_tmp$id)),
    y_obs = d_tmp$dur,
    P = ncol(X),
    X = X,
    id = d_tmp$id_idx,
    pri_s_u = 1
  )
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  foutname <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"), "sim12")
  
  # # mcmc
  f_1_1 <- m1$sample(
    ld_1, iter_warmup = 500, iter_sampling = 500,
    parallel_chains = 1, chains = 1, refresh = 10, show_exceptions = T,
    max_treedepth = 11,
    output_dir = output_dir_mcmc,
    output_basename = foutname
  )
  f_1_1$summary(variables = c("shape", "b_0", "b", "u_a", "u_r"))
  
  f_1_2 <- m1$sample(
    ld_2, iter_warmup = 500, iter_sampling = 500,
    parallel_chains = 1, chains = 1, refresh = 10, show_exceptions = T,
    max_treedepth = 11,
    output_dir = output_dir_mcmc,
    output_basename = foutname
  )
  f_1_2$summary(variables = c("shape", "b_0", "b", "u_a", "u_r"))
  # c(l_spec$shape_he, l_spec$shape_eh)
  # c(l_spec$mu_exacerb, l_spec$mu_recov)
  # c(l_spec$b_trt)
  # l_spec$b_ppfev_recov
  
  # laplace approx
  f_2_1_optim <- m1$optimize(data = ld_1, jacobian = TRUE)
  f_2_1 <- m1$laplace(data = ld_1, mode = f_2_1_optim, draws = 2000)
  f_2_1$summary(variables = c("shape", "b_0", "b", "u_a", "u_r"))
  c(l_spec$shape_eh, l_spec$mu_recov, l_spec$b_ppfev_recov, l_spec$b_trt[-1], l_spec$sd_eh)
  
  f_2_2_optim <- m1$optimize(data = ld_2, jacobian = TRUE)
  f_2_2 <- m1$laplace(data = ld_2, mode = f_2_2_optim, draws = 2000)
  f_2_2$summary(variables = c("shape", "b_0", "b", "u_a", "u_r"))
  c(l_spec$shape_eh, l_spec$mu_recov, l_spec$b_ppfev_recov, l_spec$b_trt[-1], l_spec$sd_eh)
  
  # variational inference
  f_3 <- m1$pathfinder(
    ld_1, num_paths=20,
    single_path_draws=200,
    history_size=50, max_lbfgs_iters=100,
    refresh = 0, draws = 2000)
  f_3$summary(variables = c("shape", "b_0", "b", "u_a", "u_r"))
  
  
  
  
  
  f_3 <- flexsurvreg(Surv(y, rep(1, length(y)))~trt_obs,  dist="weibullPH")
  f_3
  c(shape_eh, exp(b_0), b_trt)
  
}


get_demo_spec <- function(){
  
  
  l_spec <- list()
  
  l_spec$desc <- "ppFEV1 equal in all"
  l_spec$nsim <- 100
  l_spec$mc_cores <- 40
  
  # example trials
  l_spec$nex <- 5
  # enrolment
  l_spec$N_pt <- 600
  l_spec$pt_per_day <- 0.57
  l_spec$ramp_up_days <- 120


  
  # day of enrolment
  l_spec$t0 <- rep(1, l_spec$N_pt)
  
  
  l_spec$followup <- 365
  l_spec$rmst_eh_horizon <- 25
  
  # probability of missingness
  l_spec$pr_ymis <- 0.1
  
  # take log of mean before use in rlnorm
  l_spec$age_mean <-35
  # use as is:
  l_spec$age_sd <- 0.4
  l_spec$age_min <- 10
  l_spec$age_max <- 60
  
  # frailty parameters to link recurrences
  l_spec$sd_he <- 0.3
  l_spec$sd_eh <- 0.3 
  l_spec$rho_frailty <- -0.4
  
  l_spec$ppfev_ref <- 77.5
  l_spec$ppfev_increment <- 10
  
  
  l_spec$g_a <- 10
  l_spec$g_r <- 10
  
  # linear predictor exacerbation
  l_spec$mu_exacerb <- -16.1 # originally  -4.5
  # applied to (ppfev_baseline - reference value) / increment
  # so that at zero the linear predictor relates to the reference value 
  # and a unit change corresponds relates to a 10% increment in ppfev
  # so here the log-hazard decreases for every 10% increment in ppfev
  l_spec$b_ppfev_exacerb <- -0.2
  l_spec$shape_he <- 2.9 # originally 1.1
  
  # linear predictor recovery
  l_spec$mu_recov <- -6.5  # originally -0.5
  # as above but the log-hazard increases for every 10% increment in ppfev
  # i.e. more instantaneous risk of recovery
  l_spec$b_ppfev_recov <- 0.1
  # trt is c("soc","delay","defer")
  l_spec$b_trt <- c(0, -0.3, -0.2)
  l_spec$shape_eh <- 2.75  # originally 0.9 
  
  l_spec$trt_lab <- c("soc","delay","defer")
  l_spec$trt_active <- rep(TRUE, 3)
  names(l_spec$trt_active) <- l_spec$trt_lab
  
  l_spec$is <- 1 
  l_spec$ie <- sum(l_spec$N_pt)
  l_spec
  
}


get_scale_exacerb <- function(ppfev) {
  log_scale <- mu_exacerb + beta_ppfev_exacerb * (ppfev - ppfev_ref)
  exp(log_scale)
}


get_scale_recov <- function(ppfev = 55) {
  log_scale <- mu_recov + beta_ppfev_recov * (ppfev - ppfev_ref) + beta_trt
  exp(log_scale)
}


# Manual parameter calibration and checks for simulating exacerbation and recovery.
calibrate_weibull_ph <- function(){
  # WeibullPH
  # f(x) = amx^{a-1} exp(-m x^a)
  # F(x) = 1 - exp(-m x^a)
  # a = shape, m = scale
  # covariates included through a linear model on the log scale parameter
  # scale = exp(\mu + \beta*(ppfev - ppfev_ref) + u_i)
  # 100 <= ppfev <= 55 (approx)
  # mean is scale^(-1/shape) * Gamma(1 + 1/shape)
  
  # shape_he = 1.6
  # scale_he = 0.001
  
  shape_he = 1.1
  scale_he = 0.01
  
  cat("Days at median and upper\n")
  flexsurv::qweibullPH(p = 0.5, shape = shape_he, scale = scale_he)
  flexsurv::qweibullPH(p = 0.95, shape = shape_he, scale = scale_he)
  scale_he^(-1/shape_he) * gamma(1 + 1/shape_he)
  # sanity check
  integrand_he <- function(x, shape, scale){
    flexsurv::dweibullPH(x, shape, scale) * x
  }
  integrate(integrand_he, lower = 0, upper = Inf, shape = shape_he, scale = scale_he)
  hist(flexsurv::rweibullPH(1e5, shape = shape_he, scale = scale_he))
  plot(0:365, flexsurv::pweibullPH(0:365, shape = shape_he, scale = scale_he))
  
  # linear predictor and ppfev ref
  mu_exacerb <- -4.5
  beta_ppfev_exacerb <- -0.02
  ppfev_ref <- 77.5 
  
  
  
  get_scale_exacerb(55)
  get_scale_exacerb(100)
  # plot(1:365, pweibullPH(1:365, shape = shape_he, scale = get_scale_exacerb(66)))
  
  
  
  shape_eh = 0.7
  scale_eh = 0.14
  
  cat("Days at median and upper\n")
  flexsurv::qweibullPH(p = 0.5, shape = shape_eh, scale = scale_eh)
  flexsurv::qweibullPH(p = 0.95, shape = shape_eh, scale = scale_eh)
  # mean
  scale_eh^(-1/shape_eh) * gamma(1 + 1/shape_eh)
  hist(flexsurv::rweibullPH(1e5, shape = shape_eh))
  
  plot(0:30, flexsurv::pweibullPH(0:30, shape = shape_eh, scale = scale_eh))
  
  
  # linear predictor and ppfev ref
  mu_recov <- -1.2
  beta_ppfev_recov <- 0.02
  ppfev_ref <- 77.5 
  beta_trt <- -0.3
  
  
  
  get_scale_recov(55)
  get_scale_recov(100)
  plot(0:30, flexsurv::pweibullPH(0:30, shape = shape_eh, scale = get_scale_recov(55)))
  plot(0:30, flexsurv::pweibullPH(0:30, shape = shape_eh, scale = get_scale_recov(100)))
  
  # to these we would then need to introduce correlated frailty terms.
  
  
  
  
  
  
  
  
}

#


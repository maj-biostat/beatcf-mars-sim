
if(exists("prefix_r")){
  source(paste0(prefix_r, "/libs.R"))
  source(paste0(prefix_r, "/data.R"))
} else {
  source("./R/libs.R")
  source("./R/data.R")
}

# # re-randomisation approach
# 
# 
# if(exists("prefix_r")){
#   source(paste0(prefix_r, "/libs.R"))
#   source(paste0(prefix_r, "/data.R"))
#   source(paste0(prefix_r, "/util.R"))
# } else {
#   source("./R/libs.R")
#   source("./R/data.R")
#   source("./R/util.R")
# }
get_sim05_trial_data <- function(
    l_spec
){
  # cohort
  if(is.null(l_spec$ic)){ ic <- 1 } else { ic <- l_spec$ic }
  if(is.null(l_spec$is)){ is <- 1 } else { is <- l_spec$is }
  if(is.null(l_spec$ie)){ ie <- 1 } else { ie <- l_spec$ie }
  if(is.null(l_spec$t0)){ t0 <- 1 } else { t0 <- l_spec$t0 }
  
  trt_opts <- 1:length(l_spec$p_trt_alloc)
  
  N_ic <- length(is:ie)
  
  d_i <- data.table(
    # interim id
    ic = ic, 
    # unit id
    id = 1:N_ic,
    t0 = t0
  )
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    l_spec$age_lwr, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig, lower.tail = T)
  p_gt_age_upr <- plnorm(
    l_spec$age_upr, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(N_ic, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  baseline_age <- qlnorm(u, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig)
  d_i[, age0 := baseline_age]
  trt_alloc <- trt_opts[rep(trt_opts[as.logical(l_spec$p_trt_alloc)], len = N_ic)]
  d_i[, trt := trt_alloc]
  
  d_i[, mu0 := l_spec$b_0 + l_spec$b_la * log(age0) ]
  # what we would actually observe and therefore model
  d_i[, y0 := rnorm(.N, mu0, l_spec$sigma_mu0)]
  
  d <- CJ(
    # unit id
    id = 1:N_ic,
    # spirometry observation point - excludes baseline
    t_id = l_spec$t_sprty_obs
  )
  d <- d[d_i, on = "id"]
  
  # Build basis
  # Should generalise this but haven't got around to it. 
  # Trying to make the basis such that it is easy to enter the treatment
  # effects. At the moment it works for observation at every month 1:12 inclusive
  # with four segments. So, for example, we have a linear slope that covers month
  # 1 to 3 inclusive, then from 3 to 6, 6 to 9 and 9 to 12. 
  n_pw_seg <- nrow(l_spec$basis_ix_seg)
  k <- 1
  for (k in 1:n_pw_seg) {
    # within-interval linear ramp scaled to [0,1]
    start <- l_spec$basis_ix_seg[k, 1]
    end   <- l_spec$basis_ix_seg[k, 2]
    seg_length <- end - start
    
    d[, paste0("B", k) :=
        fifelse(t_id < start, 0,
                fifelse(t_id >= end, 1,
                        (t_id - start) / seg_length))]
  }
  
  
  #--- Compute linear predictor (mean FEV1) ---
  d[, T2 := fifelse(trt == 2, 1, 0)]
  d[, T3 := fifelse(trt == 3, 1, 0)]
  d[, mu :=
      mu0 +
      (l_spec$b_t[1] * B1 + l_spec$b_t[2] * B2 + l_spec$b_t[3] * B3 + l_spec$b_t[4] * B4) +
      (l_spec$b_T[1] * T2 + l_spec$b_T[2] * T3) +
      (l_spec$b_tT[1,1] * B1*T2 + l_spec$b_tT[1,2] * B2*T2 + l_spec$b_tT[1,3] * B3*T2 + l_spec$b_tT[1,4] * B4*T2) +
      (l_spec$b_tT[2,1] * B1*T3 + l_spec$b_tT[2,2] * B2*T3 + l_spec$b_tT[2,3] * B3*T3 + l_spec$b_tT[2,4] * B4*T3)
  ]
  
  
  
  
  # Covariance matrix for AR(1) structure only applies to the follow up
  Sigma <- outer(
    1:(length(l_spec$t_sprty_obs)), 
    1:(length(l_spec$t_sprty_obs)), 
    function(i, j){
      l_spec$sigma^2 * l_spec$rho^abs(i - j)
    } )
  
  # outcome
  # i <- 600
  for(i in 1:N_ic){
    
    d[id == i, y := as.numeric(rmvnorm(n = 1, mean = d[id == i, mu], sigma = Sigma))]
  }
  
  d[, ia := NA_integer_]
  d[, t_anlys := NA_real_]
  d[, id := rep(is:ie, each = l_spec$n_sprty_obs)]
  d[, t_fu := t0 + (t_id/l_spec$n_sprty_obs)*365]
  
  setcolorder(
    d,
    c("ic", "id", "t_id", "t0", "t_fu", "age0", 
      "trt", "mu0", "y0"))
  
  d
}




get_sim05_prototype_cfg <- function(){
  
  f_cfg_test <- file.path("./etc/sim05/cfg-sim05-v01.yml")
  cfg_test <- config::get(file = f_cfg_test)
  
  l_spec <- list()
  
  l_spec$desc <- "Test"
  l_spec$nsim <- 1
  l_spec$nex <- 1
  
  # N by analysis
  l_spec$N <- sum(cfg_test$N_pt)
  
  l_spec$pt_per_day <-   cfg_test$pt_per_day
  l_spec$ramp_up_days <-  cfg_test$ramp_up_days
  
  # number of spirometry followup (excl baseline)
  l_spec$n_sprty_obs <- cfg_test$n_sprty_obs
  l_spec$t_sprty_obs <- 1:l_spec$n_sprty_obs
  # assume linearity within intervals (over integer follow up visits)
  # basis segments - each row indexs the start/stop followup visit
  l_spec$basis_ix_seg <- matrix(
    unlist(cfg_test$basis_ix_seg),
    ncol = 2, byrow = T)
  
  
  l_spec$age_mu  <- cfg_test$age_mu
  l_spec$age_sig <- cfg_test$age_sig
  l_spec$age_lwr <- cfg_test$age_lwr
  l_spec$age_upr <- cfg_test$age_upr
  
  # trt alloc
  l_spec$p_trt_alloc <- unlist(cfg_test$trt)/length(unlist(cfg_test$trt))
  
  # parameters
  l_spec$b_0 <- cfg_test$b_0
  l_spec$b_la <- cfg_test$b_la
  l_spec$b_t <- unlist(cfg_test$b_t)
  l_spec$b_t <- sapply(1:length(l_spec$b_t), function(i) eval(parse(text = l_spec$b_t[i])))
  l_spec$b_T <- unlist(cfg_test$b_T)
  l_spec$b_tT <- matrix(
    unlist(cfg_test$b_tT),
    ncol = length(l_spec$b_t), byrow = T
  )
  
  # baseline measure
  l_spec$sigma_mu0 <- cfg_test$sigma_mu0
  
  # resid
  l_spec$sigma <- cfg_test$sigma
  l_spec$rho <- cfg_test$rho
  
  
  
  lambda = l_spec$pt_per_day
  # ramp up over x months 
  rho = function(t) pmin(t/l_spec$ramp_up_days, 1)
  
  # lambda = 0.57
  # # ramp up over x months 
  # rho = function(t) pmin(t/120, 1)
  # 
  # rr <- unlist(pblapply(1:100, cl = 4, FUN=function(ii){
  #   ttt <- get_enrol_time(sum(l_spec$N), lambda, rho)
  #   max(ttt)
  # }))
  # mean(rr) / 365
  
  # day of enrolment
  loc_t0 <- get_enrol_time(sum(l_spec$N), lambda, rho)
  
  
  l_spec$ic <- 1 # interim number
  
  # next chunk of data on pts.
  if(l_spec$ic == 1){
    # starting pt index in data
    l_spec$is <- 1
    l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
  } else {
    l_spec$is <- l_spec$ie + 1
    l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
  }
  
  # study date and follow up date
  l_spec$t0 <- loc_t0[l_spec$is:l_spec$ie]
  
  
  l_spec
}


# temporarily abandoned on 2026-02-06 with the implementation below being
# incomplete. the initial goal was to use a joint model. increasingly the
# large delay between enrol and occurrence of pex (the earliest point at 
# which the treatment can have an effect) seems problematic in the sense 
# of what is the actual treatment effect we are considering and the fact
# that the assignment doesn't align with real clinical practice.
# therefore have moved on to look into a rerand design where rand occurs at 
# the time of pex, the follow up is 14 days but we also use the long term 
# follow up as well, i.e. ppfev relative to baseline.
prototype_1_data <- function(){
  # 
  # l_spec <- get_sim05_prototype_cfg()
  # 
  # str(l_spec)
  
  N = 300
  # one year follow up
  Tmax = 1
  # time increment is 1 day
  dt <- 1/365
  
  d_pt <- data.table(
    id = 1:N
  )
  
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    6, meanlog = 2.5, sdlog = 0.75, lower.tail = T)
  p_gt_age_upr <- plnorm(
    75, meanlog = 2.5, sdlog = 0.75, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(N, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  baseline_age <- qlnorm(u, meanlog = 2.5, sdlog = 0.75)
  d_pt[, age0 := baseline_age]
  
  # equal allocation ratios for all groups (should possibly use dunnett)
  d_pt[, trt := sample(0:2, N, replace = TRUE)]
  
  
  # baseline ppfev is based on log age
  # most have pretty high ppfev1 at baseline but some are as low as ~60 ish
  d_pt[, mu0 := 120 - 14 * log(age0) ]
  summary(d_pt$mu0)
  
  beta_0     <- 90
  beta_age   <- -0.25
  beta_t     <- -2.0
  
  # main effect and time by treatment interaction
  beta_A_lvl <- c(0, -1.0, -2.0)
  beta_A_slp <- c(0, -0.5, -1.0)
  
  # instantaneous shock of pex
  alpha_PEx  <- -8
  
  # residual
  sigma_eps  <- 4
  
  
  # random intercept/slope (on time)
  Sigma_b <- matrix(c(25, -3,
                      -3, 1), 2)
  
  # state model transition intensity
  # 
  DiagrammeR("
    graph LR
      Healthy-->Exacerbation
      Exacerbation-->Death
      Healthy-->Death
             ")
  
  lambda12_0 <- 2.5    # healthy -> exacerbation
  lambda13_0 <- 0.03   # healthy -> death
  
  lambda21_0 <- 2.3    # exacerbation -> healthy 
  lambda23_0 <- 0.03   # PEx -> death
  
  
  # quantify strength of association with lung function (ties the longitudinal
  # to the state model)
  # the reason these are so low stems from the fact that we exponentiate the
  # values and because they are daily rates, i.e. the intensity on a unit 
  # interval of 1 day
  gamma12 <- -0.01
  gamma21 <- 0.01
  
  gamma13 <- -0.01; 
  gamma23 <- -0.01
  
  # effect of trt on state model transitions (beyond lung function effects)
  delta12 <- c(0, 0.02, 0.04)
  delta21 <- c(0, -0.02, -0.04)
  delta13 <- c(0, 0.01, 0.02)
  delta23 <- c(0, 0.015, 0.03)
  
  # patient-level heterogeneity
  b <- rmvnorm(N, c(0, 0), Sigma_b)
  d_pt[, `:=`(b0 = b[,1], b1 = b[,2])]
  

  # multistate simulation based on mean function for longitudinal model
  # model daily over the year 
  # multistate simulation based on mean function for longitudinal model
  # model daily over the year 
  d_state_hist <- CJ(
    id = 1:N,
    day = 0:365
  )
  d_state_hist[, state := NA_real_]
  tmp_state <- d_state_hist$state
  
  i <- j <- 1
  for (i in 1:N) {
    
    day_start <- 1
    # everyone starts off in healthy state
    cur_state <- 1
    
    # everyone starts in healthy state
    tmp_state[j] <- 1
    j <- j + 1
    
    while (day_start <= 365) {
      
      # probably need to account for state because treatment cannot
      # influence ppfev until at least one pex has occurred.
      mu_t <- d_pt$mu0[i] + 
        beta_t * (day_start/365) +
        beta_A_lvl[d_pt$trt[i] + 1] +
        # convert 365 days to range (0,1]
        beta_A_slp[d_pt$trt[i] + 1] * (day_start/365) +
        d_pt$b0[i] +
        d_pt$b1[i] * (day_start/365) +
        ifelse(cur_state == 2, alpha_PEx, 0)
      
      h12 <- lambda12_0 * exp(
        # association structure - simplest possible is:
        # current value association structure - log hazard of event at time t 
        # is linearly associated with the value of the longitudinal submodel’s 
        # linear predictor also evaluated at time t
        # for 1->2 - higher value of ppfev leads to lower transition rate
        gamma12 * mu_t +
          # individual specific covariates, just stick to direct effect of trt
          delta12[d_pt$trt[i] + 1])
      
      # for 2->1 - higher value of ppfev leads to higher transition rate
      h21 <- lambda21_0 * exp(gamma21 * mu_t +
                                delta21[d_pt$trt[i] + 1])
      
      h13 <- lambda13_0 * exp(gamma13 * mu_t +
                                delta13[d_pt$trt[i] + 1])
      
      h23 <- lambda23_0 * exp(gamma23 * mu_t +
                                delta23[d_pt$trt[i] + 1])
      
      # round(c(h12, h21, h13, h23), 3)
      
      if (cur_state == 1) {
        # pex
        p12 <- 1 - exp(-h12 * dt)
        # death
        p13 <- 1 - exp(-h13 * dt)
        # remain healthy
        p11 <- 1 - (p12 + p13)
        
        # round(c(p12, p13, p11), 3)
        
        next_state <- sample(1:3, size = 1, prob = c(p11, p12, p13))
        
      } else if (cur_state == 2) {
        # recover
        p21 <- 1 - exp(-h21 * dt)
        # death
        p23 <- 1 - exp(-h23 * dt)
        # remain in pex
        p22 <- 1 - (p21 + p23)
        
        # round(c(p21, p22, p23), 3)
        
        next_state <- sample(1:3, size = 1, prob = c(p21, p22, p23))
        
      }
      
      # state 3 is absorbing P(x = 3 | x = 3) = 1
      
      tmp_state[j] <- next_state
      j <- j + 1
      
      # update local tracking of state
      cur_state <- next_state
      day_start <- day_start + 1
      
    }
  }
  
  d_state_hist[, state := tmp_state]
  
  # now gather these daily states into start/end periods
  # episode number amounts to a cumulative sum of the times that the time 
  # ordered values for state changes by person
  d_state <- d_state_hist[
    , episode := cumsum(state != data.table::shift(state, fill = first(state))),
    by = id
  ][
    # summarise start and end time by id and episode
    , .(
      state = first(state),
      t_start = first(day),
      t_end   = last(day)
    ),
    by = .(id, episode)
  ][
    , episode := NULL
  ]
  d_state[id == 1]
  
  # on average, how many people have exacerbations
  d_state[state == 2, .N, keyby = id][
    , .(mu = mean(N), min = min(N), max = max(N))]
  
  # collate all the transitions
  d_transition <- d_state[
    , .(
      from = state,
      to   = data.table::shift(state, type = "lead")
    ),
    by = id
  ][!is.na(to)]
  
  # summarise the total number of transitions from and to each of the permitted
  # states
  d_transition[
    , .N,
    by = .(from, to)
  ][order(from, to)]
  

  # Observed Longitudinal ppfev1
  visit_times <- seq(0, 1, by = 1/12)
  
  d_visit <- CJ(
    id   = d_pt$id,
    time = visit_times
  )
  
  # convert visit time to day index
  d_visit[, day := floor(time * 365)]
  
  setkey(d_state, id, t_start, t_end)
  
  d_visit <- d_state[
    d_visit,
    on = .(
      id,
      t_start <= day,
      t_end   >= day
    ),
    nomatch = 0
  ]
  d_visit[, `:=`(t_start = NULL, t_end = NULL)]
  
  d_long <- merge(
    d_visit,
    d_pt,
    by = "id",
    all.x = TRUE
  )
  
  d_long[, mu := 
           mu0 +
           beta_t * time +
           beta_A_lvl[trt + 1] +
           beta_A_slp[trt + 1] * time +
           b0 +
           b1 * time +
           ifelse(state == 2, alpha_PEx, 0)
  ]
  d_long[, ppfev1 := rnorm(.N, mean = mu, sd = sigma_eps)]
  
  d_long[state == 2, summary(mu)]
  d_long[state == 1, summary(mu)]
  
  
  ggplot(d_long, aes(x = time, y = ppfev1, group = id)) +
    geom_point(aes(col = as.factor(state))) +
    geom_line()
  #
  
  
  list(
    d_pt = d_pt,
    d_long = d_long,
    d_state = d_state
  )
  
  
}




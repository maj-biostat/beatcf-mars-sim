


weibullPH_summary_stats <- function(shape_ph, scale_ph, tau = 25){
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
  
  
  # plot(0:365, SweibullPH(0:365, shape_ph, scale_ph))
  
  rmst <- integrate(SweibullPH, 0, tau, shape_ph = shape_ph, scale_ph = scale_ph)
  
  
  c(median = w_med, mean = w_mu, sd = w_sd, rmst = rmst$value)
}

lower_incomplete_gamma <- function(a, z) {
  pgamma(z, shape = a) * gamma(a)
}

rmst_weibull_scale <- function(tau, w_shape, w_scale) {
  s <- 1 / w_shape
  z <- (tau / w_scale)^w_shape
  (w_scale / w_shape) * gamma(s) * pgamma(z, shape = s)
}

rmst_weibull_ph <- function(tau, w_shape, w_scale) {
  s <- 1 / w_shape
  z <- w_scale * tau^w_shape
  (1 / w_shape) * w_scale^(-1 / w_shape) * gamma(s) * pgamma(z, shape = s)
}

SweibullPH <- function(x, shape_ph, scale_ph){
  1-flexsurv::pweibullPH(x, shape_ph, scale_ph)
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
  # expected value
  scale_he^(-1/shape_he) * gamma(1 + 1/shape_he)
  # sanity check for expected value
  integrand_he <- function(x, shape, scale){
    flexsurv::dweibullPH(x, shape, scale) * x
  }
  integrate(integrand_he, lower = 0, upper = Inf, shape = shape_he, scale = scale_he)
  hist(flexsurv::rweibullPH(1e5, shape = shape_he, scale = scale_he))
  plot(0:365, flexsurv::pweibullPH(0:365, shape = shape_he, scale = scale_he), type = "l")
  
  # linear predictor and ppfev ref
  mu_exacerb <- -4.5
  beta_ppfev_exacerb <- -0.2
  ppfev_ref <- 77.5 
  
  
  
  # get_scale_exacerb(55)
  # get_scale_exacerb(100)
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
  
  
  
  # get_scale_recov(55)
  # get_scale_recov(100)
  plot(0:30, flexsurv::pweibullPH(0:30, shape = shape_eh, scale = get_scale_recov(55)))
  plot(0:30, flexsurv::pweibullPH(0:30, shape = shape_eh, scale = get_scale_recov(100)))
  
  # to these we would then need to introduce correlated frailty terms.
  
  
  
  
  
  
  
  
}

#

sim_weibullPH_rmst <- function(){
  
  l_spec <- get_demo_spec()
  
  # intended to compute the rmst by treatment group marginalising over the 
  # covariate space on which the model is based.
  
  # simulate from population and then average out irrelevant terms to get 
  # rmst across trt groups
  
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    l_spec$age_min, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd, lower.tail = T)
  p_gt_age_upr <- plnorm(
    l_spec$age_max, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(1e4, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  age <- qlnorm(u, meanlog = log(l_spec$age_mean), sdlog = l_spec$age_sd)
  
  ppfev_baseline <- (ppfev_0(age, sd_ppfev = 3) - l_spec$ppfev_ref) / l_spec$ppfev_increment
  
  # indep gamma frailty but with same param values for each transition
  u_eh <- rgamma(1e4, shape = l_spec$g_a, rate = l_spec$g_r)
  
  arms = 1:3
  names(arms) <- c("soc", "defer", "discont")
  
  b_trt <- c(0, -0.25805, -0.2585)
  
  trt <- 1
  rmst_mu <- unlist(lapply(arms, function(trt) {
    
    # and then compute the recovery time on that basis
    lambda <- u_eh * exp(l_spec$mu_recov + 
                           l_spec$b_ppfev_recov * ppfev_baseline + 
                           b_trt[trt])
    
    rmst <- numeric(length(lambda))
    
    for(i in seq_along(rmst)){
      res <- integrate(
        SweibullPH, 0, 
        l_spec$rmst_eh_horizon, 
        shape_ph = l_spec$shape_eh, 
        scale_ph = lambda[i])
      
      rmst[i] <- res$value
    }
    
    mean(rmst)
    
  }))
  
  rmst_mu
  rmst_mu[2] - rmst_mu[1]
  rmst_mu[3] - rmst_mu[1]
  
}

#' Exponential decay for ppfev conditional on age with normal noise
#' 
#' f(x) = fev_max * exp(-k (age - age_ref)^p) + \epsilon
#' \epsilon \sim N(0, sd_fev)
#' 
ppfev_0 <- function(age, age_min = 10, ppfev0_max = 100, k = 0.1, p = 0.5, sd_ppfev = 0) {
  ppfev0_max * exp(-k * (age-age_min)^p) + rnorm(length(age), 0, sd_ppfev)
}



#' Inhomogeneous Poisson process
#' 
#' Simulate fixed n by thinning
sim_ipp_thinning <- function(N = 2500, lambda = 1.52,
                           rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}



#' Inhomogeneous Poisson process (Ogata’s algorithm)
#' 
#' Generate event times up to specified number of events configured to 
#' emulate an enrolment process where there is a linear ramp up to a maximum
#' rate parameter after which the process functions as a homogenous PP with 
#' fixed rate.
#' 
#' IPP has same properties as poisson process, except for the fact that its rate 
#' is a function of time, i.e., \lambda=\lambda(t)
#' 
#' For example, a football team may score goals at a higher rate at the end of 
#' a game than at the beginning if the opposing team tires more quickly. 
#' Similarly, a clinical trial may have slow recruitment at the beginning and
#' gradually ramp up to a stable rate (generally this is a bit artificial but
#' is probably better than just assuming a constant enrolment rate from the 
#' start). 
#' 
#' Define Inhomogeneous Poisson process as a counting process: {N(t), t >= 0}, 
#' so that N has integer values that never decrease over time, but jump up at 
#' random times, and we specify the following 4 conditions:
#' 
#' N(0) = 0
#' increments are independent (Markov) but not stationary
#' P(N(t + h) − N(t) = 1) = λ(t)h + o(h)
#' P(N(t + h) − N(t) > 1) = o(h)
#'
#' The number of arrivals in any interval is a Poisson random variable but 
#' the parameter can depend on the location of the interval.
#' 
#' N(t+s)−N(t) \sim \text{Poisson}(\int_t^{t+s} \lambda(a) da)
#' 
#' Ogata algorithm involves:
#' 
#' Find a constant \lambda_{max} >= \lambda(t)
#' Simulate a HPP with rate \lambda_{max}
#' Keep each event at time t with with probability \lambda(t)/\lambda_{max}
#' 
#' where \lambda_{max} is the rate at which the enrolment stabilises.
sim_ipp_ogata <- function(
    N,
    lambda = function(t, lambda_inf = 1.52, ramp_up = 90) {
      lambda_inf * pmin(t/ramp_up, 1)
    },
    lambda_max,
    ramp_up_period
) {
  t <- 0
  events <- numeric(0)
  while (length(events) < N) {
    # Propose next event from homogeneous PP
    t <- t + rexp(1, rate = lambda_max)
    
    # Accept with probability lambda(t) / lambda_max
    lambda_t <- lambda(t, lambda_max, ramp_up_period)
    if (runif(1) < lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  # event times
  events
}



#' Utility function for generating a boilerplate simulation spec used to 
#' control DGP, decision thresholds and basically anything needed in the 
#' prototyping.
#' 
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
  # trt is c("soc","defer","discont")
  l_spec$b_trt <- c(0, -0.3, -0.2)
  l_spec$shape_eh <- 2.75  # originally 0.9 
  
  
  l_spec$pri_s_u <- 0.1
  
  l_spec$trt_lab <- c("soc","defer","discont")
  l_spec$trt_active <- rep(TRUE, 3)
  names(l_spec$trt_active) <- l_spec$trt_lab
  
  l_spec$is <- 1 
  l_spec$ie <- sum(l_spec$N_pt)
  l_spec
  
}



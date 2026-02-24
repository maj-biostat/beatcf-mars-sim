

library(data.table)
library(pbapply)
library(survival)
library(alda)
library(cmdstanr)



m1 <- cmdstanr::cmdstan_model("stan/sim06-v01.stan")

simulate_ipp_fixed_n <- function(
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


# 
prototype_06_data <- function(){
  
  N <- 1000
  
  # one year follow up
  Tmax = 1
  # time increment is 1 day
  dt <- 1/365
  
  d_pt <- data.table(
    id = 1:N
  )
  
  # enrolment
  t_0 <- simulate_ipp_fixed_n(
    N = N,
    lambda = function(t, lambda_inf, ramp_up){lambda_inf * pmin(t/ramp_up, 1)},
    lambda_max = 1.52,
    ramp_up = 90
  )
  
  d_pt[, t_0 := t_0]
  
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    6, meanlog = 2.5, sdlog = 0.75, lower.tail = T)
  p_gt_age_upr <- plnorm(
    75, meanlog = 2.5, sdlog = 0.75, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(N, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  baseline_age <- qlnorm(u, meanlog = 2.5, sdlog = 0.75)
  d_pt[, age0 := baseline_age]
  
  # baseline ppfev is based on log age
  # most have pretty high ppfev1 at baseline but some are as low as ~60 ish
  d_pt[, mu0 := 120 - 14 * log(age0) ]
  
  # here we create the occurrences of pex for each pt. each pt may have a 
  # number of occurrences. we consider the duration of each episode in
  # subsequent sections of code
  
  # the occurrence of pex (not duration) is assumed to be a hpp conditional 
  # on age. this gives us the inter arrival times for pex and thus the 
  # number of pex per pt
  lambda_0 <- 0.005
  lambda_i <- lambda_0 * exp(0.05 * d_pt$mu0/100) 
  fu_dur <- 14
  l_pex <- list()
  for(i in 1:N){
    
    # create arbitrary number of event interarrival times
    t_pex <- rexp(100, lambda_i[i])
    
    # all those less than fu_dur are dropped
    t_pex <- t_pex[t_pex > fu_dur]
    # accumulate the times
    t_pex <- cumsum(t_pex)
    # stop follow up at one year
    t_pex <- t_pex[t_pex < 365]
    
    # skip those without exacerbations
    if(length(t_pex) == 0) next
    
    l_pex[[i]] <- cbind(id = d_pt$id[i], t_pex)
    
  }
  d <- data.table(do.call(rbind, l_pex))
  d[, id := as.integer(id)]
  
  d <- d_pt[d, on = .(id)]
  d[, t_pex := t_0 + t_pex]
  d[, id_pex := 1:.N, keyby = id]
  
  
  
  # at each t_pex, the pt is assigned to a treatment option
  # the pex has a duration and impact on ppfev1 conditional on trt assignment
  # for now just focus on the duration of the pex to 14 days which will be 
  # transformed into days free from symptoms to 14 days with zero also 
  d[, trt := sample(1:3, .N, replace = TRUE)]
  
  b_0 <- 2.2
  u_i <-  rnorm(length(unique(d$id)), 0, 0.07)
  d[, ix := .GRP, by = id]
  b_trt <- c(0, 0.5, -0.5)
  
  rho_ij <- exp(b_0 + u_i[d$ix] + b_trt[d$trt])
  d[, y_tru := rexp(.N, 1/rho_ij)]
  
  d[, y := pmin(y_tru, 14)]
  d[, evt := as.integer(y_tru <= 14)]
  d[, y_int := ceiling(y)]
  
  d[, ix := NULL]
  
  # table(d$y)
  d[]
  
}



run_proto_sim <- function(){
  
  
  e = NULL
  r <- pbapply::pblapply(
    X=1:10, cl = 4, FUN=function(ix) {
      
      ll <- tryCatch({
        
        
        # create data
        d <- prototype_06_data()
        
        
        
        # risk table
        f_0 <- survfit(Surv(y_int, evt)~1, data = d)
        d_rt <- data.table(
          t = f_0$time, 
          n_risk = f_0$n.risk,
          n_evt = f_0$n.event
        )
        d_rt[, n_risk_1 := data.table::shift(n_risk, fill = 0, type = "lead")]
        d_rt[, n_cens := n_risk - n_evt - n_risk_1]
        # pt at the start of the period that recovered by end
        # in discrete surv hazard is the conditional probability of 
        # event given the evt hasn't yet occurred
        d_rt[, h := n_evt / n_risk]
        d_rt[, S := f_0$surv]
        # can add first row with S = 1 and h = NA 
        d_rt
        #
        
        
        # for each row in d, expand out the days to recovery for each exacerbation
        # with an indicator as to whether there was censoring or not
        d_long <- d[rep(seq_len(.N), y_int), .(id, trt, id_pex, evt, y_int)]
        
        # crease an indicator variable for period - we need a term in the linear
        # predictor for each day
        d_long[, period := 1:.N, keyby = .(id, id_pex)]
        d_long[, evt_2 := 0]
        d_long[evt == 1, evt_2 := replace(evt_2, .N, 1), by = .(id, id_pex)]
        d_long[, evt := NULL]
        setnames(d_long, "evt_2", "evt")
        
        # sort out ids
        d_lu <- data.table(id = unique(d_long$id))
        d_lu[, id_new := 1:.N]
        
        d_long <- d_lu[d_long, on = .(id)]
        
        # 
        ld <- list(
          N = nrow(d_long),
          N_pt = nrow(d_lu),
        
          P = length(unique(d_long$period)), 
          period = d_long$period,
          evt = d_long$evt,
          
          K = length(unique(d_long$trt)),
          trt = d_long$trt,
          
          id = d_long$id_new
        )

        
        output_dir_mcmc <- paste0(getwd(), "/tmp")
        foutname <- paste0(
          format(Sys.time(), format = "%Y%m%d%H%M%S"))
        
        f_1 <- m1$sample(
          ld, iter_warmup = 1000, iter_sampling = 1000,
          parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = F,
          max_treedepth = 10,
          output_dir = output_dir_mcmc,
          output_basename = foutname
        )
        
        f_1$summary(variables = c("alpha", "s_u", "b_trt"))
        
        
        
        d_test <- copy(d_long)
        d_test[, trt := factor(trt)]
        d_test[, period := factor(period, levels = 1:14)]
        
        
        f_2 <- glm(evt ~ period + trt, 
                   family = "binomial", 
                   data = d_test)
        summary(f_2)
        
        f_3 <- inla(evt ~ -1 + period + trt + f(id_new, model = "iid"),
             family = "binomial",
             data = d_test,
             control.predictor = list(compute = TRUE)
        )
        summary(f_3)
        
        f_3$summary.fixed
        # plot(f_3$marginals.fixed$period1[, 1], f_3$marginals.fixed$period1[, 2], type = "l")
        inla.hpdmarginal(0.95, f_3$marginals.fixed$period1)
        
        
        f_2$marginals.fixed
        
        # todo - implement the above in inla
        
        # non-constant hazard - increasing hazard
        
        # rw prior on alpha
        
        # address episode effects
        
        # Past treatment can affect current health state
        # 
        # Current health state affects recovery
        
        
      },
      error=function(e) {
        message(" ERROR in MCLAPPLY LOOP " , e);
        message(traceback())
        stop(paste0("Stopping with error ", e))
      })
      
      # ll$d_all[, sum(N)]
      ll
    })
  
  
  
}
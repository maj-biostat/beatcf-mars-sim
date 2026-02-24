

library(data.table)
library(pbapply)
library(survival)
library(alda)
library(cmdstanr)
library(ggplot2)



# assume a simple multi state model where the pt oscillates between 
# healthy and exacerbation (ignore death)
simulate_patient_daily_state <- function(
    horizon = 365,
    
    # Infection / recovery parameters
    
    # daily log-odds of pex
    alpha_01 = -5.5, 
    # ppfev cov (lower ppfev1 => high risk of pex)
    beta_01  = -0.2,   
    
    # baseline recovery log-odds
    alpha_10 = -4.5,      
    # higher ppFEV1 -> faster recovery
    beta_10  =  0.1,     
    # daily increment in log-odds of recov
    gamma_10 = 0.3  ,      
    
    
    trt_effect = c(0, 0.2, -0.2),  # treatment log-odds effects (1,2,3)
    
    # ppFEV1 parameters
    ppfev_mean = 0,
    ppfev_sd   = 1,
    # both on the natural ppfev1 scale
    acute_drop = 0.1, 
    daily_recovery  = 0.02
) {
  
  # Storage
  state <- integer(horizon)
  trt   <- integer(horizon)
  ppfev <- numeric(horizon)
  p01 <- numeric(horizon)
  p10 <- numeric(horizon)
  
  # Baseline
  state[1] <- 0
  trt[1]   <- 0
  ppfev[1] <- rnorm(1, ppfev_mean, ppfev_sd)
  
  episode_duration <- 0
  
  # for this patient
  t <- 2
  for (t in 2:horizon) {
    
    if (state[t - 1] == 0) {
      
      # Healthy -> Exacerbation
      p01[t] <- inv_logit(alpha_01 + beta_01 * ppfev[t - 1])
      new_state <- rbinom(1, 1, p01[t])
      
      if (new_state == 1) {
        # pex occured, randomise treatment at onset:
        trt[t] <- sample(1:3, 1)
        
        # assume there is some acute ppFEV1 drop
        ppfev[t] <- ppfev[t - 1] + acute_drop + rnorm(1, 0, 0.01)
        
        episode_duration <- 0
        
      } else {
        
        # continue in healthy state, no trt rand
        trt[t]   <- 0
        
        # gradual (linear) decline in ppfev (in stdevs)
        ppfev[t] <- ppfev[t - 1] - 0.01/365 + rnorm(1, 0, 0.01)
      }
      
      state[t] <- new_state
      
    } else {
      # Exacerbation -> Recovery
      z <- trt[t - 1]
      p10[t] <- inv_logit(
        alpha_10 + beta_10 * ppfev[t - 1] + trt_effect[z] + gamma_10 * episode_duration
      )
      
      recovered <- rbinom(1, 1, p10[t])
      
      if (recovered == 1) {
        state[t] <- 0
        trt[t]   <- 0
        ppfev[t] <- ppfev[t - 1] + rnorm(1, 0, 0.02)
      } else {
        episode_duration <- episode_duration + 1
        state[t] <- 1
        trt[t]   <- z
        ppfev[t] <- ppfev[t - 1] + daily_recovery + rnorm(1, 0, 0.02)
      }
    }
  }
  
  d <- data.table(
    day   = 1:horizon,
    state = state,
    trt   = trt,
    ppfev = ppfev,
    p01 = p01,
    p10 = p10
  )
  
  
  d[]
}








inv_logit <- function(x) 1 / (1 + exp(-x))



n = 600

simulate_cohort_daily_state <- function(n = 600) {
  
  l <- vector("list", n)
  
  for (i in seq_len(n)) {
    d <- simulate_patient_daily_state()
    d$id <- i
    
    l[[i]] <- d
  }
  
  d_sim <- rbindlist(l)
  
  # run-length–encoding to move to start stop format:
  d_w <- d_sim[, run := rleid(state)][, .(
    start = min(day),
    stop  = max(day),
    state = state[1L]
  ),
  by = .(id, run)]
  d_w[, dur := stop - start]
  d_w[id == 76 & state == 1]
  d_w[state == 1, .N, keyby = id][, .(mu = mean(N), min = min(N), max = max(N))]
  
  
  # remove all those who do not have pex
  id_no_pex <- d_sim[, .N, keyby = .(id, state)][state == 0 & N == 365, id]
  d_sim <- d_sim[!(id %in% id_no_pex)]
  
  # remove all those with zero day duration pex
  id_zero_dur_pex <- d_w[state == 1 & dur == 0, id]
  d_sim <- d_sim[!(id %in% id_zero_dur_pex)]
  # 
  # plot(d$ppfev, main= "ppfev")
  # plot(d$state, main= "state")
  # plot(d$p01, main= "p01")
  # plot(d$p10, main= "p10")
  
  length(unique(d_sim$id))
  d_sim[]
  
  
  # ggplot(d_sim, aes(x = day, y = ppfev, group = id)) +
  #   geom_line(lwd = 0.3) +
  #   scale_y_continuous(limit = c(-4, 4))
  

  
  
}


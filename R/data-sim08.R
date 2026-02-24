

library(data.table)
library(pbapply)
library(survival)
library(alda)
library(cmdstanr)
library(ggplot2)
library(truncnorm)


# assume a simple multi state model where the pt oscillates between 
# healthy and exacerbation (ignore death)
simulate_patient_daily_state <- function(
) {
  
  N = 600
  fu_study <- 365
  fu_exac <- 14
  
  # assume continuation model for periods of health and exacerbation
  # here, the hazards are conditional probabilities.
  # \lambda(t | x) = Pr(T = t | T>= t, x) 
  # unconditional prob of failure at t is:
  # P(T = t | x) = \lambda(t|x) \prod_{s=1}^{t-1}(1 - \lambda(s|x))
  # ie the product of the conditional probabilities of not having an 
  # exacerbation in the first t-1 intervals and the probability of 
  # having an exacerbation at time t.
  
  # want a monotonically increasing hazard (log continuation ratios) for
  # the occurrence of an exacerbation in interval x (1, 2, 3, ...)
  alpha_01 <- function(x) -6.5 + 0.005 * x + 0.00001 * x^2
  
  # plot(1:fu_study, alpha_01(1:fu_study))
  
  # vectorised version of the unconditional probability for all days to x
  # compute hazard for day 1, 2, 3, ... x
  lambda_01_t <- plogis(alpha_01(1:fu_study))
  # plot(1:fu_study, lambda_01_t)
  
  # compute survival FOR ALL 0, 1, 2, ... (x-1) with surv at 0 being 1 by 
  # definition
  # For time 1: 1 (no prior times)
  # For time 2: survival_t[1] = 1 * (1 - lambda_01_t[1])
  # For time 3: survival_t[1] * 1 * (1 - lambda_01_t[1]) * (1 - lambda_01_t[2])
  # etc.
  # Can be obtained from shifted cumulative product:
  cum_survival <- cumprod(c(1, 1 - lambda_01_t[1:(fu_study-1)]))
  p_uncondit <- lambda_01_t * cum_survival
  
  # plot(1:365, p_uncondit)
  # plot(1:365, cumsum(p_uncondit))
  
  # if underflow becomes a problem due to the products then could maybe use
  # a exp log sum approach per:
  # log_lambda_t <- log(lambda_t)
  # log_cum_surv <- c(0, cumsum(log(1 - lambda_t[1:(364)])))
  # exp(log_lambda_t + log_cum_surv)
  
  # plot(1:365, p_uncondit)
  
  
  
  # want a monotonically increasing hazard (log continuation ratios) for
  # the recovery from exacerbation in interval x (1, 2, 3, ...)
  alpha_10 <- function(x) -3.5 + 0.025 * x^2
  # plot(1:fu_exac, plogis(alpha_10(1:fu_exac)))
  
  b_trt <- c(0, -1, 1)
  
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
  d_pt[, age := baseline_age]
  
  # d_pt[, ppfev_0 := 120 - 14 * log(age_0) ]
  
  d <- d_pt[CJ(id = d_pt$id, day = 1:365), on = "id"]
  d[, age := age + ((day - 1)/365)]
  d[day == 1, state := 0]
  
  i <- 1
  for(i in 1:N){
  
    state <- numeric(fu_study)
    ppfev <- numeric(fu_study)
    trt <- numeric(fu_study)
    
    state[1] <- 0
    ppfev[1] <- d[id == i & day == 1, 120 - 14 * log(age)]
    trt[1] <- 0
    
    healthy_run <- 1
    pex_run <- 1
    
    # conditional prob of pex 
    j <- 2
    for(j in 2:fu_study){
      
      if(state[j-1] == 0){
        
        pex_run <- 1
        healthy_run <- healthy_run + 1
        # healthy -> pex (or remain in healthy)
        p_01 <- plogis(alpha_01(healthy_run) + 0.01 * ppfev[j-1])
        
        # draw pex evt in interval j given healthy up to and including j-1
        state[j] <- rbinom(1, 1, p_01)
        
        # update ppfev - gradually declining over the year
        ppfev[j] <- d[id == i & day == j, 120 - 14 * log(age) + rnorm(1, 0, 0.2)]
        trt[j] <- 0
        
      } else if(state[j-1] == 1){
        
        pex_run <- pex_run + 1
        healthy_run <- 1
        
        # randomise trt
        trt[j] <- sample(0:2, size = 1)
        
        # pex -> healthy (or remain in pex) 
        p_10 <- plogis(alpha_10(j) + 0.01 * ppfev[j-1] + b_trt[trt[j]])
       
        # draw health transition in interval j given pex up to and including j-1
        state[j] <- rbinom(1, 1, p_10)
        
        # update ppfev
        
        
        
      }
      
      
      
      
      
      
      
      
      
    }
    
    
    
    # conditional probabilities for this individual
    lambda_t <- plogis(alpha_01(1:fu_study) + 0.2 * d_pt[i, ppfev_0])
    
    tt <- 1
    while(tt < 365){
      
      # draw event at time tt given survival to tt
      pex <- rbinom(1, 1, lambda_t[tt])
      
      
      
      
    }
    
    
  }
  
  
  
  d[]
}








inv_logit <- function(x) 1 / (1 + exp(-x))




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




# source dependencies
toks <- unlist(data.table::tstrsplit(getwd(), "/")) 
if(toks[length(toks)] == "beatcf-mars-sim"){
  prefix_cfg <- "./etc/sim16/"
  prefix_stan <- "./stan"
  prefix_fig <- "./fig"
  prefix_data <- "./data"
  prefix_r <- "./R"
} else {
  prefix_cfg <- "../etc/sim16/"
  prefix_stan <- "../stan"
  prefix_fig <- "../fig"
  prefix_data <- "../data"
  prefix_r <- "../R"
}

source(paste0(prefix_r, '/libs.R'))
source(paste0(prefix_r, '/init.R'))
source(paste0(prefix_r, '/util.R'))



#' 
get_sim18_pt <- function(
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
  
  eh_u <- rnorm(1, 0, l_spec$u_sd_eh)
  he_u <- rnorm(1, 0, l_spec$u_sd_he)
  
  # make time to first exacerbation exponential
  
  
  # duration of exacerbation on treatment is log-normal
  
  
  # time to recrudescence for a proportion of the cohort
  
  
  
  # entry time is zero, in a healthy state
  day_of_fu <- 0
  state <- "H"
  v_state <- rep("H", l_spec$decision_fu)
  
  exac_count  <- 0
  
  names(l_spec$b_trt_eh) <- l_spec$trt_lab
  names(l_spec$b_trt_he) <- l_spec$trt_lab
  
  # days to first exacerbatin (in healthy state)
  days_in_init <- 0
  t_exac <- 0
  # evt indicator
  status <- 0
  
  # trt isn't assigned unless an exac occurs
  trt <- "none"
  
  repeat {
    
    days_in_init <- days_in_init + 1
    bin_ix <- l_spec$v_lu_he_bin[days_in_init]
    lp <- l_spec$a_he[bin_ix] +
      he_u +
      l_spec$b_ppfev_he * ppfev_0
    
    event <- rbinom(1, 1, plogis(lp))
    
    if(event == 1){
      t_exac = days_in_init
      status <- 1
      state <- "E"
      exac_count <- 1
      break
    }
    
    if(days_in_init >= 365){
      t_exac <- l_spec$followup
      status <- 0
      break
    }
  }
  c(t_exac, status, state)
  
  if(status == 1){
    
    trt <- sample(l_spec$trt_lab[l_spec$trt_active], 1)
    days_in_state <- 0
    
    for(i in 1:l_spec$decision_fu){
      
      if (state == "E") {
        
        bin_ix <- l_spec$v_lu_eh_bin[days_in_state + 1L]
        lp <- l_spec$a_eh[bin_ix] + eh_u +
          l_spec$b_ppfev_eh * ppfev_0 +
          l_spec$b_trt_eh[trt]
        
        days_in_state <- days_in_state + 1
        y <- rbinom(1,1,plogis(lp))
        
        if(y == 1){
          state <- "H"
          days_in_state <- 0
        }
        
      } else {
        
        bin_ix <- l_spec$v_lu_he_bin[days_in_state + 1L]
        # carryover trt is already set
        lp <- l_spec$a_he[bin_ix] + he_u +
          l_spec$b_ppfev_he * ppfev_0 +
          l_spec$b_trt_he[trt]
        
        days_in_state <- days_in_state + 1
        y <- rbinom(1,1,plogis(lp))
        
        if(y == 1){
          state <- "E"
          exac_count <- exac_count + 1
          days_in_state <- 0
        }
      }
      
      v_state[i] <- state
    }
    
  }
  
  # the last observation should probably be censored because we do not see
  # the completion of the exacerbation or healthy spell
  d_pt <- data.table(
    day_of_exac = t_exac,
    evt = status,
    trt = trt,
    exac_days = sum(v_state == "E"),
    exac_count = exac_count
  )
  d_pt <- cbind(age = age, ppfev_0 = ppfev_0, d_pt)
  d_pt[, age := age + day_of_fu/l_spec$followup]
  d_pt[]
}

#' Wrapper to invoke state transition simulation for individual patients used
#' to create cohort
#' 
get_sim18_cohort <- function(l_spec){
  
  id_cohort <- l_spec$is:l_spec$ie
  N_cohort <- length(id_cohort)
  
  d_cohort <- data.table(
    id  = l_spec$is:l_spec$ie,
    trt = sample(l_spec$trt_lab[l_spec$trt_active], N_cohort, replace = T)
  )
  
  # Draw correlated standard normal scores (copula scale)
  # (w1, w2) ~ BVN(0, 0, 1, 1, rho)
  w_1 <- rnorm(N_cohort)
  # indep:
  z_2 <- rnorm(N_cohort)
  # incorporate correlation
  w_2 <- l_spec$rho * w_1 + sqrt(1 - l_spec$rho^2) * z_2
  # should be approx l_spec$rho
  # cor(w_1, w_2)
  
  # Map margins to uniform via standard normal CDF
  u_1 <- pnorm(w_1)
  u_2 <- pnorm(w_2)
  
  # Map uniforms through the inverse log-logistic AFT CDF 
  # Log-logistic AFT: 
  # the logit transform converts the [0,1] variable into logistic
  # log(T) = mu + sigma * logit(u),  
  # u ~ Uniform(0,1)

  # The cdf for T is the same as the cdf for epsi where epsi = logit(u) 
  # P(T <= t)= P(epsi <= z) = u = F(t)
  # logit(u) = (log(t)-mu)/sigma
  # log(t) = mu + sigma*logit(u)
  
  mu_1 <- l_spec$b_1_0 + l_spec$b_1_trt[d_cohort[, trt]]
  mu_2 <- l_spec$b_2_0 + l_spec$b_2_trt[d_cohort[, trt]]
  
  logit <- function(p) log(p / (1 - p))
  
  true_t1 <- exp(mu_1 + l_spec$sig_1 * qlogis(u_1))
  # gap time AFTER recovery, i.e.
  true_t2 <- exp(mu_2 + l_spec$sig_2 * qlogis(u_2))   
  
  # calendar time of next exacerbation is true_t1 + true_t2
  
  d_cohort[, `:=`(true_t1 = true_t1, true_t2 = true_t2)]
  
  # Administrative censoring on the calendar scale:
  # Pt followed from randomisation (assume to be time 0) up to l_spec$followup days
  # T1 is censored if true_t1 > l_spec$followup days. If T1 is observed,
  # T2's clock starts at true_t1 (calendar time), and the remaining budget
  # for observing T2 is (l_spec$followup days - true_t1); T2 is censored if
  # true_t2 exceeds that remaining budget. 
  
  # If T1 itself is censored, T2 is
  # structurally unobserved/undefined (the pt never recovered within
  # the window - unlikely).
  # In that situation T2 as right-censored at a small positive placeholder 
  # gap time (1e-6 days; see t2_placeholder below). 
  
  # This is not an arbitrary hack: for the Gaussian-copula
  # log-logistic likelihood, the both-censored contribution P(T1>t1, T2>t2) 
  # collapses to S1(t1) alone as t2 -> 0+ (since F2(t2) -> 0 and the copula CDF
  # C(u1,u2;rho) -> 0 along with it), so this placeholder correctly
  # encodes "T2's clock never started" as "pure T1 censoring" without
  # needing a separate code path in the Stan model. 
  
  
  d_cohort[, evt_1 := as.integer(true_t1 <= l_spec$followup)]
  d_cohort[, t1_obs := pmin(true_t1, l_spec$followup)]
  
  d_cohort[, remaining_window := l_spec$followup - t1_obs]
  d_cohort[, evt_2 := as.integer(evt_1 == 1L & true_t2 <= remaining_window)]
  
  ## Placeholder gap time for T2 when T1 itself is censored: this needs to
  ## be small relative to the realistic T2 scale (so the both-censored
  ## likelihood contribution correctly collapses to ~S1(t1), see the
  ## analytic check in the accompanying notes) WITHOUT being so close to
  ## machine epsilon that log(t2) becomes an extreme value feeding into the
  ## Stan likelihood's log-logistic transform -- 1e-6 days is many orders
  ## of magnitude below any realistic gap time in this trial while numerically
  ## benign (log(1e-6) ~= -13.8, well within normal double precision range).
  t2_placeholder <- 1e-6
  
  d_cohort[, t2_obs := fifelse(
    evt_1 == 1L,
    pmin(true_t2, pmax(remaining_window, t2_placeholder)),
    t2_placeholder   ## d2 is forced to 0 in this branch by construction
  )]
  
  d_cohort[]
  
}



#' Convert sample data.table into lists suitable for stan models
#' 
get_sim18_stan_data <- function(dd, l_spec){
  
  
  dd[trt == "def", defer := 1]
  dd[trt != "def", defer := 0]
  dd[trt == "dis", discont := 1]
  dd[trt != "dis", discont := 0]
  
  X <- model.matrix(~-1 + ppfev_0 + defer + discont, data = dd)
  
  ld <- list(
    N = dd[, .N],
    
    P = ncol(X),
    X = X,
    
    y = dd$exac_days,
    
    trt_defer_col = 2,
    trt_discont_col = 3,
    
    # zero centred effects with sd
    pri_b_0 = l_spec$pri_b_0,
    pri_b = l_spec$pri_b,
    pri_s = l_spec$pri_s,
    
    prior_only = 0
  )
  
  ld
}

#' Convert the day wise data into segments 
#' 
sim18_long_to_wide <- function(dd){
  
  dd_w <- dd[, .(
    tstart = min(day_of_fu)-1,
    tstop = max(day_of_fu),
    evt = max(y),
    state = state[1],
    trt = trt[1],
    bin = bin[1],
    ppfev_0 = ppfev_0[1]
  ), keyby = .(id, rlgrp)]
  
  # View(d_cohort_w)
  
  dd_w[, .N, keyby = .(id, state)][, .(mu = mean(N)), keyby = state]
  
  dd_w[, state := factor(state, levels = c("H", "E"))]
  
  dd_w[trt == "def", defer := 1]
  dd_w[trt != "def", defer := 0]
  dd_w[trt == "dis", discont := 1]
  dd_w[trt != "dis", discont := 0]
  
  dd_w[tstart == 0, i_entry := 1]
  dd_w[tstart != 0, i_entry := 0]
  
  dd_w[, len_seg := tstop - tstart]
  
  dd_w
}






update_sim18_cfg <- function(l_spec){
  
  
  l_spec$pri_b_0 <- unlist(l_spec$pri_b_0)
  l_spec$pri_b <- unlist(l_spec$pri_b)
  l_spec$pri_s <- unlist(l_spec$pri_s) 
  
  
  # initially all trt arms are active
  l_spec$trt_lab <- unlist(l_spec$trt_lab)
  # *** Has to be converted to logical otherwise you are just going to be 
  # indexing trt 1
  l_spec$trt_active <- as.logical(l_spec$trt_active)
  names(l_spec$trt_active) <- l_spec$trt_lab
  
  
  l_spec$par_names <- c(
    "b_1_0", 
    "b_1_trt[2]", "b_1_trt[3]",
    "b_2_0", 
    "b_2_trt[2]", "b_2_trt[3]"
  )
  
  l_spec$b_1_0 <- log(l_spec$e_b_1_0) 
  l_spec$b_1_trt <- log(unlist(l_spec$e_b_1_trt))
  names(l_spec$b_1_trt) <- l_spec$trt_lab
  
  l_spec$b_2_0 <- log(l_spec$e_b_2_0) 
  l_spec$b_2_trt <- log(unlist(l_spec$e_b_2_trt))
  names(l_spec$b_2_trt) <- l_spec$trt_lab
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$nsim)
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$nsim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
  
  l_spec
}






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



# todo 
# review betancourt warping paper etc


#' State transitions for patient to followup time
#' 
get_sim16_pt <- function(
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
  ppfev_0 <- (f_ppfev_0(age, sd_ppfev = 3) - l_spec$ppfev_ref) / l_spec$ppfev_increment
  
  eh_u <- rnorm(1, 0, l_spec$u_sd_eh)
  he_u <- rnorm(1, 0, l_spec$u_sd_he)
  
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
get_sim16_cohort <- function(l_spec){
  pt_list <- list()
  id_cohort <- l_spec$is:l_spec$ie
  i <- 1
  for(i in seq_along(id_cohort)){
    pt_list[[i]] <- cbind(id = id_cohort[i],  t0 = l_spec$t0[i], get_sim16_pt(l_spec))
  }
  
  d_cohort <- rbindlist(pt_list)
  
  d_cohort[]
}



#' Convert sample data.table into lists suitable for stan models
#' 
get_sim16_stan_data <- function(dd, l_spec){
  
  
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
sim16_long_to_wide <- function(dd){
  
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






cfg_update <- function(l_spec){
  
  # recovery bins (up to 25 days)
  l_spec$eh_bins <- unlist(l_spec$eh_bins) 
  l_spec$he_bins <- unlist(l_spec$he_bins) 
  
  l_spec$a_he <- unlist(l_spec$a_he) 
  l_spec$a_eh <- unlist(l_spec$a_eh) 
  
  
  l_spec$pri_b_0 <- unlist(l_spec$pri_b_0)
  l_spec$pri_b <- unlist(l_spec$pri_b)
  l_spec$pri_s <- unlist(l_spec$pri_s) 
  
  # trt alloc - balanced over number of trts
  
  # exacerbation -> healthy - linpred for exacerbation state impacting duration of recovery
  l_spec$b_trt_eh <- unlist(l_spec$b_trt_eh)
  l_spec$n_trt_eh <- length(l_spec$b_trt_eh)
  
  # healthy -> exacerbation - linpred for healthy state impacting period to relapse occurs
  l_spec$b_trt_he <- unlist(l_spec$b_trt_he)
  l_spec$n_trt_he <- length(l_spec$b_trt_he)
  
  names(l_spec$b_trt_eh) <- l_spec$trt_lab
  names(l_spec$b_trt_he) <- l_spec$trt_lab
  
  l_spec$par_names_pre <- c("b_0", "b")
  l_spec$par_names <- c(
    "b_0", 
    "b_ppfev0", 
    paste0("b_", l_spec$trt_lab[2:3])
  )
  
  # initially all trt arms are active
  l_spec$trt_lab <- unlist(l_spec$trt_lab)
  
  
  # *** Has to be converted to logical otherwise you are just going to be 
  # indexing trt 1
  l_spec$trt_active <- as.logical(l_spec$trt_active)
  names(l_spec$trt_active) <- l_spec$trt_lab
  
  
  # bin lookup - avoid findInterval
  l_spec$d_lu_he_bin <- data.table(
    day = l_spec$he_bins, ix_bin = seq_along(l_spec$he_bins))
  d_grid <- data.table(day = 0:max(l_spec$he_bins))
  l_spec$d_lu_he_bin <- l_spec$d_lu_he_bin[d_grid, on = "day", roll = T]
  # essential to set key otherwise this will be painfully slow
  setkey(l_spec$d_lu_he_bin, day)
  
  # optimisation 
  # vector based lookup
  l_spec$v_lu_he_bin <- l_spec$d_lu_he_bin$ix_bin
  # could just get to this as the bin boundary + 1...
  l_spec$rle_he <- rle(l_spec$v_lu_he_bin)
  l_spec$he_starts <- cumsum(c(1L, head(l_spec$rle_he$lengths, -1)))
  
  l_spec$d_lu_eh_bin <- data.table(
    day = l_spec$eh_bins, ix_bin = seq_along(l_spec$eh_bins))
  d_grid <- data.table(day = 0:max(l_spec$eh_bins))
  l_spec$d_lu_eh_bin <- l_spec$d_lu_eh_bin[d_grid, on = "day", roll = T]
  # essential to set key otherwise this will be painfully slow
  setkey(l_spec$d_lu_eh_bin, day)
  
  l_spec$v_lu_eh_bin <- l_spec$d_lu_eh_bin$ix_bin
  l_spec$rle_eh <- rle(l_spec$v_lu_eh_bin)
  l_spec$eh_starts <- cumsum(c(1L, head(l_spec$rle_eh$lengths, -1)))
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$nsim)
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$nsim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
  
  l_spec
}




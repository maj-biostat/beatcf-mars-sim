

# source dependencies
toks <- unlist(data.table::tstrsplit(getwd(), "/")) 
if(toks[length(toks)] == "beatcf-mars-sim"){
  prefix_cfg <- "./etc/sim15/"
  prefix_stan <- "./stan"
  prefix_fig <- "./fig"
  prefix_data <- "./data"
  prefix_r <- "./R"
} else {
  prefix_cfg <- "../etc/sim15/"
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
get_sim15_pt <- function(
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
  v_day_of_fu <- integer(l_spec$followup)
  v_state <- integer(l_spec$followup)
  v_y <- integer(l_spec$followup)
  v_day_in_state <- integer(l_spec$followup)
  v_trt <- character(l_spec$followup)
  ix_state <- 0
  
  exac_count  <- 0
  
  names(l_spec$b_trt_eh) <- l_spec$trt_lab
  names(l_spec$b_trt_he) <- l_spec$trt_lab
  
  while (day_of_fu < l_spec$followup) {
    
    days_in_state <- 0
    
    repeat{
      
      ix_state <- ix_state + 1
      
      if (state == "H") {
        bin_ix <- l_spec$v_lu_he_bin[days_in_state + 1L] 
        if(exac_count == 0){
          trt <- "none"
          lp <- l_spec$a_he[bin_ix] + he_u + 
            l_spec$b_ppfev_he * ppfev_0
        } else {
          # carryover trt is already set
          lp <- l_spec$a_he[bin_ix] + he_u + 
            l_spec$b_ppfev_he * ppfev_0 + 
            l_spec$b_trt_he[trt]
        }
        
      } else {
        bin_ix <- l_spec$v_lu_eh_bin[days_in_state + 1L] 
        exac_count <- exac_count + 1
        lp <- l_spec$a_eh[bin_ix] + eh_u + 
          l_spec$b_ppfev_eh * ppfev_0 + 
          l_spec$b_trt_eh[trt]
      }
      
      days_in_state <- days_in_state + 1
      
      y <- rbinom(1,1,plogis(lp))
      
      v_day_of_fu[ix_state] <- day_of_fu + days_in_state
      v_state[ix_state] <- state
      v_y[ix_state] <- y
      v_day_in_state[ix_state] <- days_in_state
      v_trt[ix_state] <- trt
      
      # c(state, days_in_state, y)
      if (y==1 || (day_of_fu+days_in_state)>=l_spec$followup) break
      
    }
    
    day_of_fu <- day_of_fu + days_in_state
    
    if (state=="H" && y==1) {
      
      state <- "E"
      # randomise treatment at the start of the exacerbation we just hit
      trt_options <- l_spec$trt_lab[l_spec$trt_active]
      trt <- sample(trt_options, 1)
    } else if (state=="E" && y==1) {
      state <- "H"
    }
    
  }
  
  # the last observation should probably be censored because we do not see
  # the completion of the exacerbation or healthy spell
  d_pt <- data.table(
    day_of_fu = v_day_of_fu,
    state = v_state,
    y = v_y,
    day_in_state = v_day_in_state,
    trt = v_trt
  )
  d_pt <- cbind(age = age, ppfev_0 = ppfev_0, d_pt)
  d_pt[, age := age + day_of_fu/l_spec$followup]
  d_pt
  
}

#' Wrapper to invoke state transition simulation for individual patients used
#' to create cohort
#' 
get_sim15_cohort <- function(l_spec){
  pt_list <- list()
  id_cohort <- l_spec$is:l_spec$ie
  i <- 1
  for(i in seq_along(id_cohort)){
    pt_list[[i]] <- cbind(id = id_cohort[i],  t0 = l_spec$t0[i], get_sim15_pt(l_spec))
  }
  
  d_cohort <- rbindlist(pt_list)
  
  d_cohort[]
}



#' Convert sample data.table into lists suitable for stan models
#' 
get_sim15_stan_data <- function(dd, l_spec){
  
  dd_w <- sim15_long_to_wide(dd)
  
  X <- model.matrix(~-1 + ppfev_0 + defer + discont, data = dd_w)
  
  ld <- list(
    N   = dd_w[, .N],
    N_id = dd_w[, .N, keyby = id][, .N],
    id   = dd_w$id,
    y = dd_w$evt,
    # 1 is healthy, 2 is exacerbation
    state = dd_w[, as.integer(state)],
    bin = dd_w$bin,
    N_he_bin = length(l_spec$a_he),
    N_eh_bin = length(l_spec$a_eh),
    
    i_entry = dd_w$i_entry,
    len_seg = dd_w$len_seg,
    
    P = ncol(X),
    X = X,
    trt_defer_col = 2,
    trt_discont_col = 3,
    
    # zero centred effects with sd
    pri_a_he_mu = l_spec$pri_a_he_mu,
    pri_a_he_s = l_spec$pri_a_he_s,
    
    pri_a_eh_mu = l_spec$pri_a_eh_mu,
    pri_a_eh_s = l_spec$pri_a_eh_s,
    
    pri_b_s = l_spec$pri_b_s,
    
    pri_u_he_r = l_spec$pri_u_he_r,
    pri_u_eh_r = l_spec$pri_u_eh_r,
    
    prior_only = 0
  )
  
  stopifnot(all(ld$bin <= max(c(ld$N_he_bin, ld$N_eh_bin))))
  # str(ld)
  
  ld = ld
}

#' Convert the day wise data into segments 
#' 
sim15_long_to_wide <- function(dd){
  
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

# expected time in exacerbation state
sim_policy <- function(
    a_he, b_he_ppfev, b_he_trt, u_he, 
    a_eh, b_eh_ppfev, b_eh_trt, u_eh, 
    ppfev_0, policy = 1, l_spec) {
  
  t <- 0L
  # start in the exacerbation state
  state <- 2L   # 1 = H, 2 = E
  # we start in an exacerbation state so 
  # this starts at 1 not zero.
  exac_count <- 1L
  
  total_days_in_E <- 0L
  
  # precompute
  lp_ppfev_he <- b_he_ppfev * ppfev_0
  lp_ppfev_eh <- b_eh_ppfev * ppfev_0
  
  decision_fu <- l_spec$decision_fu
  
  rle_he <- l_spec$rle_he
  rle_eh <- l_spec$rle_eh
  
  while (t < decision_fu) {
    
    days_in_state <- 0L
    
    # track position in RLE
    bin_pos <- 1L
    offset_in_bin <- 0L
    
    repeat {
      
      if (state == 1L) {
        
        bin_ix <- rle_he$values[bin_pos]
        bin_width <- rle_he$lengths[bin_pos] - offset_in_bin
        
        if (exac_count == 0L) {
          lp <- a_he[bin_ix] + u_he + lp_ppfev_he
        } else {
          lp <- a_he[bin_ix] + u_he + lp_ppfev_he + b_he_trt[policy]
        }
        
      } else {
        
        bin_ix <- rle_eh$values[bin_pos]
        bin_width <- rle_eh$lengths[bin_pos] - offset_in_bin
        
        lp <- a_eh[bin_ix] + u_eh + lp_ppfev_eh + b_eh_trt[policy]
      }
      
      # convert to probability
      p <- plogis(lp)
      
      # constant hazard within piecewise bins
      # draw geometric (number of days until event)
      # rgeom gives failures before success => +1
      wait <- rgeom(1, p) + 1L
      
      if (wait <= bin_width) {
        # event happens in this bin
        days_in_state <- days_in_state + wait
        
        if (state == 2L) {
          total_days_in_E <- total_days_in_E + 
            min(wait, decision_fu - t)
        }
        
        break
        
      } else {
        # consume full bin and move to next
        days_in_state <- days_in_state + bin_width
        
        if (state == 2L) {
          total_days_in_E <- total_days_in_E + 
            min(bin_width, decision_fu - t)
        }
        
        bin_pos <- bin_pos + 1L
        offset_in_bin <- 0L
      }
      
      if ((t + days_in_state) >= decision_fu) break
    }
    
    t <- t + days_in_state
    
    if (t >= decision_fu) break
    
    # state transition
    if (state == 1L) {
      state <- 2L
      exac_count <- exac_count + 1L
    } else {
      state <- 1L
    }
  }
  
  total_days_in_E
}

# B_max = 100
# N_pt = 100
# ppfev_baseline = unique(d_cohort$ppfev_baseline)
calc_trt_effect <- function(
    d_post, B_max = 100, N_pt = 100, 
    # unique ppfev0 from sample
    ppfev_0, l_spec){
  
  B <- min(B_max, nrow(d_post))
  
  m_a_he <- as.matrix(d_post[, .SD, .SDcols = patterns("a_he")])
  v_b_he_ppfev <- d_post$b_he_1
  m_b_he_trt <- cbind(0, d_post$b_he_2, d_post$b_he_3)
  v_u_he <- rnorm(B, 0, d_post$u_sd_he)
  m_a_eh <- as.matrix(d_post[, .SD, .SDcols = patterns("a_eh")])
  v_b_eh_ppfev <- d_post$b_eh_1
  m_b_eh_trt <- cbind(0, d_post$b_eh_2, d_post$b_eh_3)
  v_u_eh <- rnorm(B, 0, d_post$u_sd_eh)
  
  v_ppfev_0 <- sample(ppfev_0, size = B, replace = T)

  m_soc <- matrix(NA, nrow = B, ncol = 2)
  m_def <- matrix(NA, nrow = B, ncol = 2)
  m_dis <- matrix(NA, nrow = B, ncol = 2)
  
  soc <- numeric(N_pt)
  def <- numeric(N_pt)
  dis <- numeric(N_pt)
  
  trt_ix <- setNames(seq_along(l_spec$trt_lab), l_spec$trt_lab)
  
  i <- 1
  for(i in 1:B){
    
    for(j in 1:N_pt){
      soc[j] <- sim_policy(
        a_he = m_a_he[i, ],
        b_he_ppfev = v_b_he_ppfev[i],
        b_he_trt = m_b_he_trt[i, ],
        u_he = v_u_he[i],
        
        a_eh = m_a_eh[i, ],
        b_eh_ppfev = v_b_eh_ppfev[i],
        b_eh_trt = m_b_eh_trt[i, ],
        u_eh = v_u_eh[i],
        
        ppfev_0 = v_ppfev_0[i],
        policy = trt_ix[1], 
        l_spec
      )
    }
    
    for(j in 1:N_pt){
      def[j] <- sim_policy(
        a_he = m_a_he[i, ],
        b_he_ppfev = v_b_he_ppfev[i],
        b_he_trt = m_b_he_trt[i, ],
        u_he = v_u_he[i],
        
        a_eh = m_a_eh[i, ],
        b_eh_ppfev = v_b_eh_ppfev[i],
        b_eh_trt = m_b_eh_trt[i, ],
        u_eh = v_u_eh[i],
        
        ppfev_0 = v_ppfev_0[i],
        policy = trt_ix[2], 
        l_spec
      )
    }
    
    for(j in 1:N_pt){
      dis[j] <- sim_policy(
        a_he = m_a_he[i, ],
        b_he_ppfev = v_b_he_ppfev[i],
        b_he_trt = m_b_he_trt[i, ],
        u_he = v_u_he[i],
        
        a_eh = m_a_eh[i, ],
        b_eh_ppfev = v_b_eh_ppfev[i],
        b_eh_trt = m_b_eh_trt[i, ],
        u_eh = v_u_eh[i],
        
        ppfev_0 = v_ppfev_0[i],
        policy = trt_ix[3], 
        l_spec
      )
    }
    
    # d_fig <- data.table(soc, def, dis)
    # d_fig <- melt(d_fig, measure.vars = names(d_fig))
    # ggplot(d_fig, aes(x = value, group = variable, col = variable)) + geom_density()
    
    # health state given by decision fu less the exacerbation time
    m_soc[i, ] = colMeans(cbind(soc, l_spec$decision_fu - soc))
    m_def[i, ] = colMeans(cbind(def, l_spec$decision_fu - def))
    m_dis[i, ] = colMeans(cbind(dis, l_spec$decision_fu - dis))
    
  }
  
  d_res <- data.table(
    mu_soc_H = m_soc[, 2], 
    mu_def_H = m_def[, 2], 
    mu_dis_H = m_dis[, 2], 
    mu_soc_E = m_soc[, 1], 
    mu_def_E = m_def[, 1], 
    mu_dis_E = m_dis[, 1])
  
  d_res[, `:=`(
    del_def_H = mu_def_H - mu_soc_H,
    del_dis_H = mu_dis_H - mu_soc_H,
    del_def_E = mu_def_E - mu_soc_E,
    del_dis_E = mu_dis_E - mu_soc_E
  )]
  
  d_res
  
}




cfg_update <- function(l_spec){
  
  # recovery bins (up to 25 days)
  l_spec$eh_bins <- unlist(l_spec$eh_bins) 
  l_spec$he_bins <- unlist(l_spec$he_bins) 
  
  l_spec$a_he <- unlist(l_spec$a_he) 
  l_spec$a_eh <- unlist(l_spec$a_eh) 
  
  # trt alloc - balanced over number of trts
  
  # exacerbation -> healthy - linpred for exacerbation state impacting duration of recovery
  l_spec$b_trt_eh <- unlist(l_spec$b_trt_eh)
  l_spec$n_trt_eh <- length(l_spec$b_trt_eh)
  
  # healthy -> exacerbation - linpred for healthy state impacting period to relapse occurs
  l_spec$b_trt_he <- unlist(l_spec$b_trt_he)
  l_spec$n_trt_he <- length(l_spec$b_trt_he)
  
  names(l_spec$b_trt_eh) <- l_spec$trt_lab
  names(l_spec$b_trt_he) <- l_spec$trt_lab
  
  l_spec$par_names_pre <- c("a_he", "b_he", "u_sd_he", "a_eh", "b_eh", "u_sd_eh")
  l_spec$par_names <- c(
    paste0("a_he_", seq_along(l_spec$a_he)),
    paste0("b_he_", seq_along(c(l_spec$b_ppfev_he, l_spec$b_trt_he[-1]))),
    "u_sd_he",
    paste0("a_eh_", seq_along(l_spec$a_eh)),
    paste0("b_eh_", seq_along(c(l_spec$b_ppfev_eh, l_spec$b_trt_eh[-1]))),
    "u_sd_eh"
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



example_sim15_v02 <- function(){
  
  ##########
  # Create dataset of arbitrary size
  
  source("R/libs.R")
  source("R/init.R")
  source("R/util.R")
  source("R/data-sim15.R")
  
  f_cfgsc <- file.path("./etc/sim15/cfg-sim15-v04.yml")
  l_spec <- config::get(file = f_cfgsc)
  
  l_spec <- cfg_update(l_spec)
  
  pt_list <- list()
  
  id_cohort <- 1:10000
  i <- 1
  for(i in seq_along(id_cohort)){
    pt_list[[i]] <- cbind(id = id_cohort[i],  t0 = NA, get_sim15_pt(l_spec))
  }
  d_cohort <- rbindlist(pt_list)
  d_cohort[, .N, keyby = trt]

  ###########
  # Distribution of exacerbations
  
  d_cohort[state == "H", bin := l_spec$v_lu_he_bin[day_in_state + 1L]]
  d_cohort[state == "E", bin := l_spec$v_lu_eh_bin[day_in_state + 1L]]
  d_cohort[, rlgrp := rleid(id, state, trt, bin)]
  d_w <- sim15_long_to_wide(d_cohort)
  
  # Identify distinct periods in each state within each patient
  # Creates running episode index w/in pt
  d_w[, per_id := cumsum(
    # identify shift in state between this and the next rec
    state != data.table::shift(
      state, fill = data.table::first(state))) + 1L, 
    keyby = .(id)
  ]
  
  d_fig <- d_w[, .(
    n_E = uniqueN(per_id[state == "E"])
  ), keyby = .(id)
  ]
  
  d_fig[, .N, keyby = .(n_E)]
  
  # Number of exacerbations
  # https://github.com/tidyverse/ggplot2/issues/2051
  ggplot(d_fig, aes(x = n_E)) +
    geom_bar(aes(y = after_stat(prop), group = 1)) +
    scale_x_continuous("Exacerbations", breaks = 0:10) +
    scale_y_continuous("Proportion of Pts", 
                       breaks = seq(0, 1, by = 0.1)) 
  
  
  ###########
  # Duration of first healthy period
  
  d_fig <- d_w[state == "H" & per_id == 1]
  d_fig <- d_fig[, .(dur_H_1 = sum(len_seg)), keyby = id]
  
  # Duration of first healthy period
  mean(d_fig$dur_H_1 == 365) # 20% are healthy throughout
  summary(d_fig$dur_H_1)  
  ggplot(d_fig, aes(x = dur_H_1)) +
    geom_histogram(bins = 30) +
    scale_x_continuous("Duration of first H episode") 
  
  
  ###########
  # Days spent in the exacerbation state
  
  d_fig <- d_w[
    state == "E", .(dur_eh = sum(len_seg)), keyby = .(id, per_id, trt)]
  
  d_trt <- d_fig[, .(mu = mean(dur_eh)), keyby = trt]
  d_trt[trt == "def", mu] - d_trt[trt == "soc", mu]
  l_spec$dec_eh_delta_ni
  
  ggplot(d_fig, aes(x = dur_eh)) +
    geom_histogram(bins = 20) +
    scale_x_continuous("Duration (days)")  +
    facet_wrap(~trt)
  
  ############
  # Repeat sim - Days spent in the exacerbation state 
  
  f_cfgsc <- file.path("./etc/sim15/cfg-sim15-v04.yml")
  l_spec <- config::get(file = f_cfgsc)
  l_spec <- cfg_update(l_spec)
  n_sim <- 400
  delta_eh_def <- as.vector(n_sim)
  id_cohort <- 1:500
  i <- j <- 1
  
  for(j in 1:n_sim){
    
    pt_list <- list()
    for(i in seq_along(id_cohort)){
      pt_list[[i]] <- cbind(id = id_cohort[i],  t0 = NA, get_sim15_pt(l_spec))
    }
    d_cohort <- rbindlist(pt_list)
    
    d_cohort[state == "H", bin := l_spec$v_lu_he_bin[day_in_state + 1L]]
    d_cohort[state == "E", bin := l_spec$v_lu_eh_bin[day_in_state + 1L]]
    d_cohort[, rlgrp := rleid(id, state, trt, bin)]
    d_w <- sim15_long_to_wide(d_cohort)
    
    d_w[, per_id := cumsum(
      # identify shift in state between this and the next rec
      state != data.table::shift(
        state, fill = data.table::first(state))) + 1L, 
      keyby = .(id)
    ]
    
    d_res <- d_w[
      state == "E", .(dur_eh = sum(len_seg)), keyby = .(id, per_id, trt)]
    d_res <- d_res[, .(mu = mean(dur_eh)), keyby = trt]
    delta_eh_def[j] <- d_res[trt == "def", mu] - d_res[trt == "soc", mu]
  }
  hist(delta_eh_def, main = sprintf("Delta mu = %.2f b_trt_eh[2] = %.2f", mean(delta_eh_def), l_spec$b_trt_eh[2]))
  abline(v = mean(delta_eh_def), col = 2)
  
  
}



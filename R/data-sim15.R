

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
      
      # states_list[[length(states_list)+1]] <- data.table(
      #   day_of_fu= day_of_fu + days_in_state,
      #   state=state,
      #   y=y,
      #   day_in_state=days_in_state,
      #   trt=trt
      # )
      
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
    i <- i + 1
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
    
    pri_sd_he = 4,
    pri_sd_eh = 4
  )
  # str(ld)
  
  ld
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
  
  dd_w[trt == "defer", defer := 1]
  dd_w[trt != "defer", defer := 0]
  dd_w[trt == "discont", discont := 1]
  dd_w[trt != "discont", discont := 0]
  
  dd_w[tstart == 0, i_entry := 1]
  dd_w[tstart != 0, i_entry := 0]
  
  dd_w[, len_seg := tstop - tstart]
  
  dd_w
}

# expected exacerbation days for a randomly selected participant from the trial population under a given treatment policy
sim_policy_1 <- function(a_he, b_he, u_he, 
                       a_eh, b_eh, u_eh, 
                       # individual ppfev0 from our sample
                       ppfev_0,
                       # policy which we want to apply over the fu
                       policy = 1, l_spec) {
  
  t <- 0
  state <- "H"
  exac_count <- 0
  
  total_days_in_H <- 0
  total_days_in_E <- 0
  
  lp_ppfev_he <- l_spec$b_ppfev_he * ppfev_0
  lp_ppfev_eh <- l_spec$b_ppfev_eh * ppfev_0
  
  while (t < l_spec$followup) {
    # days in whatever state you happen to be in
    days_in_state <- 0
    repeat{
      
      if (state == "H") {
        bin_ix <- l_spec$v_lu_he_bin[days_in_state + 1L]
        if(exac_count == 0){
          lp <- a_he[bin_ix] + u_he + 
            lp_ppfev_he
        } else {
          # carryover trt is already set
          lp <- a_he[bin_ix] + u_he + 
            lp_ppfev_he + b_he[policy]
        }
        
      } else {
        bin_ix <- l_spec$v_lu_eh_bin[days_in_state + 1L] 
        exac_count <- exac_count + 1
        lp <- a_eh[bin_ix] + u_eh + 
          lp_ppfev_eh + b_eh[policy]
      }
      
      # end of interval marks the transition, for example if you had
      # (0, 1], (1, 2], (2, 3] in a health state H and then had a transition in
      # (3, 4] to E, you would have had a total of 3 days in the healthy state
      # because the transition is assumed to occur on the last interval
      
      days_in_state <- days_in_state + 1
      
      # cat(paste(days_in_state, " days in state ", state, "\n"))
      
      if (state=="H" && (t+days_in_state) <= l_spec$followup) {
        total_days_in_H <- total_days_in_H + 1
      }
      if (state=="E" && (t+days_in_state) <= l_spec$followup) {
        total_days_in_E <- total_days_in_E + 1
      }
      
      y <- rbinom(1, 1, plogis(lp))
      
      # stop repeat if we had an event (transition to other state) or reached fu
      if (y==1 || (t+days_in_state)>=l_spec$followup) break
      
    }
    
    t <- t + days_in_state
    
    if (state=="H" && y==1) {
      state <- "E"
    } else if (state=="E" && y==1) {
      state <- "H"
    }
    
  }
  
  # c(total_days_in_H, total_days_in_E, total_days_in_H + total_days_in_E)
    
  total_days_in_E
}

sim_policy_2 <- function(a_he, b_he, u_he, 
                            a_eh, b_eh, u_eh, 
                            ppfev_0,
                            policy = 1,
                            l_spec) {
  
  t <- 0L
  state <- 1L   # 1 = H, 2 = E
  exac_count <- 0L
  
  total_days_in_E <- 0L
  
  # precompute
  lp_ppfev_he <- l_spec$b_ppfev_he * ppfev_0
  lp_ppfev_eh <- l_spec$b_ppfev_eh * ppfev_0
  
  followup <- l_spec$followup
  
  # bin lookup vectors
  v_he <- l_spec$v_lu_he_bin
  v_eh <- l_spec$v_lu_eh_bin
  
  max_bin_he <- length(v_he)
  max_bin_eh <- length(v_eh)
  
  while (t < followup) {
    
    days_in_state <- 0L
    
    repeat {
      
      # current bin index
      k <- days_in_state + 1L
      
      if (state == 1L) max_bin <- max_bin_he  
      if (state == 2L) max_bin <- max_bin_eh  
      
      if (state == 1L && k > max_bin_he) k <- max_bin  # safety
      if (state == 2L && k > max_bin_eh) k <- max_bin  # safety
      
      if (state == 1L) {
        
        bin_ix <- v_he[k]
        
        if (exac_count == 0L) {
          lp <- a_he[bin_ix] + u_he + lp_ppfev_he
        } else {
          lp <- a_he[bin_ix] + u_he + lp_ppfev_he + b_he[policy]
        }
        
      } else {
        
        bin_ix <- v_eh[k]
        
        lp <- a_eh[bin_ix] + u_eh + lp_ppfev_eh + b_eh[policy]
      }
      
      # convert to probability
      p <- 1 / (1 + exp(-lp))
      
      # draw geometric (number of days until event)
      # rgeom gives failures before success → +1
      wait <- rgeom(1, p) + 1L
      
      # how many days remain in this bin?
      # (since bins are defined implicitly via lookup)
      # we detect next bin change
      next_k <- k
      
      while (next_k <= max_bin && 
             ((state == 1L && v_he[next_k] == bin_ix) ||
              (state == 2L && v_eh[next_k] == bin_ix))) {
        next_k <- next_k + 1L
      }
      
      bin_width <- next_k - k
      
      if (wait <= bin_width) {
        # event occurs in this bin
        days_in_state <- days_in_state + wait
        
        if (state == 2L) {
          total_days_in_E <- total_days_in_E + 
            min(wait, followup - t)
        }
        
        break
        
      } else {
        # no event in this bin => jump to next bin
        days_in_state <- days_in_state + bin_width
        
        if (state == 2L) {
          total_days_in_E <- total_days_in_E + 
            min(bin_width, followup - t)
        }
      }
      
      if ((t + days_in_state) >= followup) break
    }
    
    t <- t + days_in_state
    
    if (t >= followup) break
    
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
  m_b_he <- as.matrix(d_post[, .SD, .SDcols = patterns("b_he")])
  v_u_he <- rnorm(B, 0, d_post$u_sd_he)
  m_a_eh <- as.matrix(d_post[, .SD, .SDcols = patterns("a_eh")])
  m_b_eh <- as.matrix(d_post[, .SD, .SDcols = patterns("a_eh")])
  v_u_eh <- rnorm(B, 0, d_post$u_sd_eh)
  
  v_ppfev_0 <- sample(ppfev_0, size = B, replace = T)

  v_soc <- numeric(B)
  v_def <- numeric(B)
  v_dis <- numeric(B)
  
  soc <- numeric(N_pt)
  def <- numeric(N_pt)
  dis <- numeric(N_pt)
  
  trt_ix <- setNames(seq_along(l_spec$trt_lab), l_spec$trt_lab)
  
  i <- 1
  for(i in 1:B){
    
    for(j in 1:N_pt){
      soc[j] <- sim_policy_2(
        a_he = m_a_he[i, ],
        b_he = m_b_he[i, ],
        u_he = v_u_he[i],
        
        a_eh = m_a_eh[i, ],
        b_eh = m_b_eh[i, ],
        u_eh = v_u_eh[i],
        
        ppfev_0 = v_ppfev_0[i],
        policy = trt_ix[1], 
        l_spec
      )
    }
    
    for(j in 1:N_pt){
      def[j] <- sim_policy_2(
        a_he = m_a_he[i, ],
        b_he = m_b_he[i, ],
        u_he = v_u_he[i],
        
        a_eh = m_a_eh[i, ],
        b_eh = m_b_eh[i, ],
        u_eh = v_u_eh[i],
        
        ppfev_0 = v_ppfev_0[i],
        policy = trt_ix[2], 
        l_spec
      )
    }
    
    for(j in 1:N_pt){
      dis[j] <- sim_policy_2(
        a_he = m_a_he[i, ],
        b_he = m_b_he[i, ],
        u_he = v_u_he[i],
        
        a_eh = m_a_eh[i, ],
        b_eh = m_b_eh[i, ],
        u_eh = v_u_eh[i],
        
        ppfev_0 = v_ppfev_0[i],
        policy = trt_ix[3], 
        l_spec
      )
    }
    
    
    
    v_soc[i] = mean(soc)
    v_def[i] = mean(def)
    v_dis[i] = mean(dis)
    
  }
  
  d_res <- data.table(soc = v_soc, def = v_def, dis = v_dis)
  d_res[, `:=`(
    delta_def = def - soc,
    delta_dis = dis - soc
  )]
  
  d_res
  
}






example_sim15_v02 <- function(){
  
  
  pt_list <- list()
  
  id_cohort <- 1:600
  i <- 1
  for(i in seq_along(id_cohort)){
    pt_list[[i]] <- cbind(id = id_cohort[i],  t0 = NA, get_sim15_pt(l_spec))
    i <- i + 1
  }
  
  d_cohort <- rbindlist(pt_list)
  
  d_cohort[, .N, keyby = trt]
  
  
  d_cohort[state == "H", bin := findInterval(d, l_spec$he_bins)]
  d_cohort[state == "E", bin := findInterval(d, l_spec$eh_bins)]
  # constant intervals (start stop) defined by
  d_cohort[, grp := rleid(id, state, trt, bin)]
  
  View(d_cohort)
  
  # days_E <- d_cohort[state == "E", .N, keyby = .(id, bin)][, N]
  # hist(days_E)
  
  # min(t) - 1 to get intervals respecting (start, stop] 
  d_cohort_w <- d_cohort[, .(
    tstart = min(t)-1,
    tstop = max(t),
    evt = max(y),
    state = state[1],
    trt = trt[1],
    bin = bin[1],
    ppfev0 = ppfev_baseline[1]
  ), keyby = .(id, grp)]
  
  View(d_cohort_w)
  
  d_cohort_w[, .N, keyby = .(id, state)][, .(mu = mean(N)), keyby = state]
  
  d_cohort_w[, state := factor(state, levels = c("H", "E"))]
  
  d_cohort_w[trt == "defer", defer := 1]
  d_cohort_w[trt != "defer", defer := 0]
  
  d_cohort_w[trt == "discont", discont := 1]
  d_cohort_w[trt != "discont", discont := 0]
  
  d_cohort_w[tstart == 0, i_entry := 1]
  d_cohort_w[tstart != 0, i_entry := 0]
  
  d_cohort_w[, len_seg := tstop - tstart]
  
  View(d_cohort_w)
  
  X <- model.matrix(~-1 + ppfev0 + defer + discont, data = d_cohort_w)
  
  
  ld <- list(
    N   = nrow(d_cohort_w),
    N_id = length(unique(d_cohort_w$id)),
    id   = d_cohort_w$id,
    y = d_cohort_w$evt,
    # 1 is healthy, 2 is exacerbation
    state = as.integer(d_cohort_w$state),
    bin = d_cohort_w$bin,
    N_he_bin = length(l_spec$mu_he_0),
    N_eh_bin = length(l_spec$mu_eh_0),
    
    i_entry = d_cohort_w$i_entry,
    len_seg = d_cohort_w$len_seg,
    
    P = ncol(X),
    X = X,
    trt_defer_col = 2,
    trt_discont_col = 3,
    
    pri_sd_he = 4,
    pri_sd_eh = 4
  )
  
  
  m1 <- cmdstanr::cmdstan_model("stan/sim15-v02.stan")
  
  f_1 <- m1$sample(
    ld, iter_warmup = 1000, iter_sampling = 1000,
    parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = F,
    max_treedepth = 11
  )
  f_1$summary(variables = c(
    "a_he", "a_eh", "b_he", "b_eh", "sd_he", "sd_eh"
  ))
  
  d_post <- data.table(
    f_1$draws(
      format = "matrix",
      variables = c("a_he", "b_he", "sd_he",
                    "a_eh", "b_eh", "sd_eh"))
  )
  
  
  f_2_optim <- m1$optimize(data = ld, jacobian = TRUE, refresh = 0)
  f_2 <- m1$laplace(data = ld, mode = f_2_optim, draws = 2000, refresh = 0)
  
  f_2$summary(variables = c(
    "a_he", "a_eh", "b_he", "b_eh", "sd_he", "sd_eh"
  ))
}


# test_nbinom <- function(x_max = 20, mu = 8){
#   x <- 0:x_max
#   y <- pnbinom(x, 1, mu = mu)
#   plot(x, y, ylim = c(0, 1), type = "l")
#   abline(h = 1)
# }
# test_nbinom(x_max = 365, mu = 50)

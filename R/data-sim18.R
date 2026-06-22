

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





#' Wrapper to invoke state transition simulation for individual patients used
#' to create cohort
#' 
sim18_cohort <- function(l_spec){
  
  id_cohort <- l_spec$is:l_spec$ie
  N_cohort <- length(id_cohort)
  
  # produce data for every day and then chop it down to what we observe
  d_cohort <- data.table(
    id  = l_spec$is:l_spec$ie,
    day = rep(0:l_spec$max_day, N_cohort),
    trt = sample(l_spec$trt_lab[l_spec$trt_active], N_cohort, replace = T)
  )
  setorder(d_cohort, id, day)
  d_cohort[, t_0 := l_spec$t_0[id]]
  
  d_cohort[, state := NA_integer_]
  d_cohort[day==0, state:= sample(l_spec$state_opts, .N, replace = TRUE, prob = l_spec$p_init)]
  
  for(i in seq_len(N_cohort))
  {
    
    rows <- which(d_cohort$id == id_cohort[i])
    
    for(j in 2:length(rows))
    {
      
      prev <- d_cohort$state[rows[j-1]]
      
      # all intervals are structurally 1 because we simulate daily and then 
      # subset to the survey days
      interval <- 1
      
      tt <- d_cohort$day[rows[j]]
      trt <- d_cohort$trt[rows[j]]
      
      lp <-
        # treatment
        l_spec$b_trt[trt] +
        # previous state
        l_spec$b_prev[prev] +
        # time (quadratic)
        l_spec$b_time_1*tt +
        l_spec$b_time_2*tt^2 +
        # time by treatment (linear)
        l_spec$b_trt_time[trt]*tt +
        # gap length - structurally zero and is only relevant in model
        # where we have gaps in observation
        l_spec$b_gap[interval]
      
      # not very efficient but it will do
      d_cohort$state[rows[j]] <- sim18_rord_pom_3(lp,l_spec$alpha)
      
    }
    
  }
  
  d_obs <- d_cohort[day %in% l_spec$visit_days]
  
  list(
    d_cohort = d_cohort,
    d_obs = d_obs
  )
  
}



#' Convert sample data.table into lists suitable for stan models
#' 
sim18_stan_data <- function(dd, l_spec){
  
  
  setorder(dd, id, day)
  
  dd[, `:=`(
    prev_state = data.table::shift(state, 1L),
    prev_day   = data.table::shift(day, 1L)
  ), by = id]
  
  # For the day zero of onset, we have no prev state. We can either assume that they
  # were well the day before or just drop that observation and include it in day 1
  # of follow up. I do the latter.
  dd[, gap_len := day - prev_day]
  dd <- dd[day != 0]
  # gives 16 days of follow up per pt 14 whole days and then 1 obs in wk 3 and 4
  dd[, trt_idx := match(trt, l_spec$trt_lab)]
  
  dd[, x_time := copy(day)]
  # scaling makes this sample a lot faster but is a pain in the arse
  # for linear predictors with interactions
  # additionally means you need to back scale intercepts
  dd[, x_time := scale(x_time)]
  dd[, x_gap_time := factor(fifelse(dd$gap_len == 1, 1, 2))]
  dd[, x_trt := factor(trt_idx)]
  dd[, x_prev := factor(prev_state, levels = l_spec$state_opts)]
  
  X <- model.matrix(~ x_trt + x_prev +
                      x_time + I(x_time^2) +
                      x_gap_time +
                      x_trt * x_time, 
                    data = dd)
  X_mod <- X[, -1]
  
  ld <- list(
    N  = nrow(dd),
    P = ncol(X_mod),
    X = X_mod,
    y = dd$state,
    # Indexing parameters within design matrix.
    # First list any variable with no dependency on time
    ix_trt_2 = 1,
    ix_trt_3 = 2,
    ix_prev_2 = 3,
    ix_prev_3 = 4,
    ix_time_1 = 5,
    ix_time_2 = 6,
    ix_gap = 7,
    ix_trt_time_2 = 8,
    ix_trt_time_3 = 9,
    mu_days = mean(dd$day),
    sd_days = sd(dd$day)
  )
  stopifnot(names(X_mod)[ld$ix_trt_2] == "x_trt2")
  stopifnot(names(X_mod)[ld$ix_trt_3] == "x_trt3")
  stopifnot(names(X_mod)[ld$ix_prev_2] == "x_prev2")
  stopifnot(names(X_mod)[ld$ix_prev_3] == "x_prev3")
  stopifnot(names(X_mod)[ld$ix_time_1] == "x_time")
  stopifnot(names(X_mod)[ld$ix_time_2] == "I(x_time^2)")
  stopifnot(names(X_mod)[ld$ix_gap] == "x_gap_time2")
  stopifnot(names(X_mod)[ld$ix_trt_time_2] == "x_trt2:x_time")
  stopifnot(names(X_mod)[ld$ix_trt_time_3] == "x_trt3:x_time")
  
  ld
}



# construct 3 state transition probs based on linear predictor and 
# cutpoints from ordinal model
sim18_rord_pom_3 <- function(lp, alpha)
{
  p0 <- plogis(alpha[1] - lp)
  p1 <- plogis(alpha[2] - lp) - p0
  p2 <- 1 - p0 - p1
  
  sample(1:3, 1, prob = c(p0, p1, p2)) 
}

sim18_transition_matrix <- function(day, gap_ix = 1, trt, l_spec)
{
  
  P <- matrix(0,3,3)
  
  for(prev in l_spec$state_opts){
    
    lp <-
      l_spec$b_trt[trt] +
      l_spec$b_prev[prev] +
      l_spec$b_time_1 * day +
      l_spec$b_time_2 * day^2 +
      l_spec$b_trt_time[trt] * day +
      l_spec$b_gap[gap_ix]
    
    p1 <- plogis(l_spec$alpha[1] - lp)
    p2 <- plogis(l_spec$alpha[2] - lp) - p1
    p3 <- 1 - p1 - p2
    
    P[prev,] <- c(p1,p2,p3)
    
  }
  
  # from 
  rownames(P) <- l_spec$state_lab
  # to 
  colnames(P) <- l_spec$state_lab
  P
  
}




sim18_sop <- function(days = 1:28, l_spec)
{
  
  # has to start at day 1
  stopifnot(days[1] == 1)
  
  out <- matrix(NA, length(days) + 1, 3)
  
  out[1,] <- l_spec$p_init
  
  # starting point
  pi <- l_spec$p_init
  
  d_sop <- data.table()
  
  for(trt in l_spec$trt_lab){
    
    for(i in seq_along(days)){
      
      day = days[i]
      
      if(day == 1){
        gap_ix = 1;
      } 
      
      if(day > 1){
        if(days[i] - days[i-1] == 1){
          gap_ix = 1
        } else {
          gap_ix = 2
        }
      }
      
      Pt <- sim18_transition_matrix(day, gap_ix, trt, l_spec)
      
      pi <- drop(pi %*% Pt)
      
      out[i + 1,] <- pi
      
    }
    
    d_sop <- rbind(
      d_sop, 
      data.table(
        trt = trt,
        day = c(0, days),
        none = out[,1],
        mild = out[,2],
        severe = out[,3]
      )
    )
  }
  
  d_sop
  
}


update_sim18_cfg <- function(l_spec){
  
  l_spec$trt_lab <- unlist(l_spec$trt_lab)
  l_spec$state_lab <- unlist(l_spec$state_lab)
  
  l_spec$alpha <- unlist(l_spec$alpha)
  l_spec$b_trt <- unlist(l_spec$b_trt)
  names(l_spec$b_trt) <- l_spec$trt_lab
  
  # *** Has to be converted to logical otherwise you are just going to be 
  # indexing trt 1
  l_spec$trt_active <- as.logical(l_spec$trt_active)
  names(l_spec$trt_active) <- l_spec$trt_lab
  
  l_spec$state_opts <- seq_along(l_spec$state_lab)
  names(l_spec$state_opts) <- l_spec$state_lab
  
  l_spec$visit_days <- unlist(l_spec$visit_days)
  l_spec$max_day <- max(l_spec$visit_days)
  
  l_spec$b_prev <- unlist(l_spec$b_prev)
  
  l_spec$b_trt_time <- unlist(l_spec$b_trt_time)
  names(l_spec$b_trt_time) <- l_spec$trt_lab
  
  l_spec$b_gap <- unlist(l_spec$b_gap)
  l_spec$p_init <- unlist(l_spec$p_init)
  
  l_spec$smry_pars <- c("a", "b_trt", "b_prev", "b_time_1", "b_time_2", "b_gap", "b_trt_time")
  
  l_spec$full_pars <- c("a[1]", "a[2]",
                        "b_trt[1]", "b_trt[2]", "b_trt[3]",
                        "b_prev[1]", "b_prev[2]", "b_prev[3]",
                        "b_time_1", "b_time_2",
                        "b_gap[1]", "b_gap[2]",
                        "b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")
    
    
    c("a", "b_trt", "b_prev", "b_time_1", "b_time_2", "b_gap", "b_trt_time")
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$nsim)
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$nsim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
  
  l_spec
}




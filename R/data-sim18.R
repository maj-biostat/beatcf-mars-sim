

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
sim18_cohort <- function(l_spec){
  
  id_cohort <- l_spec$is:l_spec$ie
  N_cohort  <- length(id_cohort)
  n_day     <- l_spec$max_day + 1
  
  ## subject-level quantities
  trt_subj <- sample(
    l_spec$trt_lab[l_spec$trt_active],
    N_cohort,
    replace = TRUE
  )
  
  ## cohort data
  d_cohort <- data.table(
    id  = rep(id_cohort, each = n_day),
    day = rep(0:l_spec$max_day, times = N_cohort),
    trt = rep(trt_subj, each = n_day)
  )
  
  d_cohort[, t_0 := l_spec$t_0[id]]
  
  ## work directly on vectors
  state <- integer(nrow(d_cohort))
  
  ## initial states
  state[seq(1L, by = n_day, length.out = N_cohort)] <-
    sample(
      l_spec$state_opts,
      N_cohort,
      replace = TRUE,
      prob = l_spec$p_init
    )
  
  ## cache coefficients
  b_trt        <- l_spec$b_trt
  b_prev       <- l_spec$b_prev
  b_prev_time  <- l_spec$b_prev_time
  b_trt_time <- l_spec$b_trt_time
  alpha       <- l_spec$alpha
  # gap_effect  <- l_spec$b_gap[1]
  
  ## precompute time effect
  day_vals <- 0:l_spec$max_day
  time_eff <- l_spec$b_time_1 * day_vals + l_spec$b_time_2 * day_vals^2
  i <- 1
  
  for(i in seq_len(N_cohort)) {
    
    first <- (i - 1L) * n_day + 1L
    last  <- first + n_day - 1L
    
    trt_i <- trt_subj[i]
    
    for(r in (first + 1L):last) {
      
      prev <- state[r - 1L]
      tt   <- d_cohort$day[r]
      
      lp <-
        b_trt[trt_i] +
        b_prev[prev] +
        time_eff[tt + 1L] +
        b_prev_time[prev] * tt + 
        b_trt_time[trt_i] * tt 
      
      p0 <- plogis(alpha[1] - lp)
      p1 <- plogis(alpha[2] - lp)
      
      u <- runif(1)
      
      state[r] <-
        if (u < p0) {
          1L
        } else if (u < p1) {
          2L
        } else {
          3L
        }
    }
  }
  
  d_cohort[, state := state]
  
  d_obs <- copy(d_cohort)
  setorder(d_obs, id, day)
  
  # Need to compute previous state and previous day before we subset
  # otherwise we will end up referencing day 7 as the previous day from day 14
  d_obs[, `:=`(
    prev_state = data.table::shift(state, 1L),
    prev_day   = data.table::shift(day, 1L)
  ), by = id]
  
  d_obs <- d_obs[day %in% l_spec$visit_days]
  
  list(
    d_cohort = d_cohort,
    d_obs    = d_obs
  )
  
}



#' Convert sample data.table into lists suitable for stan models
sim18_stan_data <- function(dd, l_spec){
  
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
  
  dd[, x_trt := factor(trt_idx)]
  dd[, x_prev := factor(prev_state, levels = l_spec$state_opts)]
  
  X <- model.matrix(~ x_trt + x_prev +
                      x_time + I(x_time^2) +
                      x_prev * x_time +
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
    ix_prev_time_2 = 7,
    ix_prev_time_3 = 8,
    ix_trt_time_2 = 9,
    ix_trt_time_3 = 10,
    mu_days = mean(dd$day),
    sd_days = sd(dd$day)
  )
  stopifnot(names(X_mod)[ld$ix_trt_2] == "x_trt2")
  stopifnot(names(X_mod)[ld$ix_trt_3] == "x_trt3")
  stopifnot(names(X_mod)[ld$ix_prev_2] == "x_prev2")
  stopifnot(names(X_mod)[ld$ix_prev_3] == "x_prev3")
  stopifnot(names(X_mod)[ld$ix_time_1] == "x_time")
  stopifnot(names(X_mod)[ld$ix_time_2] == "I(x_time^2)")
  stopifnot(names(X_mod)[ld$ix_prev_time_2] == "x_prev2:x_time")
  stopifnot(names(X_mod)[ld$ix_prev_time_3] == "x_prev3:x_time")
  stopifnot(names(X_mod)[ld$ix_trt_time_2] == "x_trt2:x_time")
  stopifnot(names(X_mod)[ld$ix_trt_time_3] == "x_trt3:x_time")
  
  ld
}



sim18_transition_matrix <- function(day, trt, l_spec)
{
  
  P <- matrix(0,3,3)
  
  for(prev in l_spec$state_opts){
    
    lp <-
      l_spec$b_trt[trt] +
      l_spec$b_prev[prev] +
      l_spec$b_time_1 * day +
      l_spec$b_time_2 * day^2 +
      l_spec$b_prev_time[prev] * day +
      l_spec$b_trt_time[trt] * day 
    
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
  d_sop <- data.table()
  
  for(trt in l_spec$trt_lab){
    # starting point
    pi <- l_spec$p_init
    
    for(i in seq_along(days)){
      
      day = days[i]
      Pt <- sim18_transition_matrix(day, trt, l_spec)
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
  
  l_spec$b_prev_time <- unlist(l_spec$b_prev_time)
  
  # l_spec$b_gap <- unlist(l_spec$b_gap)
  
  # l_spec$b_trt_gap <- unlist(l_spec$b_trt_gap)
  # names(l_spec$b_trt_gap) <- l_spec$trt_lab
  
  l_spec$p_init <- unlist(l_spec$p_init)
  
  l_spec$smry_pars <- c(
    "a", "b_trt", "b_prev", "b_time_1", "b_time_2", 
    # "b_gap", 
    "b_prev_time", 
    "b_trt_time")
  
  l_spec$full_pars <- c("a[1]", "a[2]",
                        "b_trt[1]", "b_trt[2]", "b_trt[3]",
                        "b_prev[1]", "b_prev[2]", "b_prev[3]",
                        "b_time_1", "b_time_2",
                        # "b_gap[1]", "b_gap[2]", "b_gap[3]",  "b_gap[4]",
                        "b_prev_time[1]", "b_prev_time[2]", "b_prev_time[3]",
                        "b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")
  
  l_spec$non_zero_pars <- c("a[1]", "a[2]",
                        "b_trt[2]", "b_trt[3]",
                        "b_prev[2]", "b_prev[3]",
                        "b_time_1", "b_time_2",
                        # "b_gap[2]", "b_gap[3]",  "b_gap[4]",
                        "b_prev_time[2]", "b_prev_time[3]",
                        "b_trt_time[2]", "b_trt_time[3]")
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$nsim)
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$nsim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
    
  # convert to the induced expected difference in sojourn times
  l_spec$delta_lab <- c("delta_def", "delta_dis")
  d_sop <- sim18_sop(1:28, l_spec)
  d_tbl <- dcast(d_sop, day ~ trt, value.var = "none")
  
  l_spec$dur_tru <- c( 
    d_tbl[day > 0, sum(soc)],
    d_tbl[day > 0, sum(def)],
    d_tbl[day > 0, sum(dis)]
  )
  names(l_spec$dur_tru) <- l_spec$trt_lab
  
  l_spec$delta_dur_tru <- c( 
    d_tbl[day > 0, sum(def - soc)],
    d_tbl[day > 0, sum(dis - soc)]
  )
  names(l_spec$delta_dur_tru) <- l_spec$delta_lab
   
  l_spec
}

# use to approximate treatment effects that will lead to desired differences
# across arms in terms of soujourn time
sim18_calibrate_trt <- function(l_spec){
  
  source("R/libs.R")
  source("R/init.R")
  source("R/util.R")

  f_cfgsc <- file.path("./etc/sim18/cfg-sim18-v01.yml")
  l_spec <- config::get(file = f_cfgsc)

  l_spec <- update_sim18_cfg(l_spec)
  l_spec$t_0 <- seq_along(1:sum(l_spec$N_pt))

  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$N_pt)
  
  # where are we at the moment
  l_spec$b_trt
  days <- 1:max(l_spec$visit_days)
  d_sop <- sim18_sop(days, l_spec)
  dcast(d_sop[day > 0], day ~ trt, value.var = "none" )[day %in% c(1, 4, 7, 14, 21, 28)]
  
  
  l_spec$dec_delta_ni <- 1
  explr_interval <- c(-5, 5)
  message("Traget NI margin  : ", l_spec$dec_delta_ni)
  
  f_obj <- function(b_trt) {
    l_spec$b_trt["def"] <- b_trt
    d_sop <- sim18_sop(days, l_spec)
    d_tbl <- dcast(d_sop[day > 0], day ~ trt, value.var = "none")
    delta <- d_tbl[, sum(def - soc)]
    (delta + l_spec$dec_delta_ni)^2
  }
  
  f_trt <- stats::optimise(f = f_obj, interval = explr_interval)
  # nominated treatment effect
  message("Trt effect to induce NI  : ", f_trt$minimum)
  
  # Recompute as a sanity check
  l_spec$b_trt["def"] <- f_trt$minimum
  d_sop <- sim18_sop(days, l_spec)
  
  d_tbl <- dcast(d_sop[day > 0], day ~ trt, value.var = "none" )
  # Duration in no symptom state by trt
  d_tbl[, .(def = sum(def), dis = sum(dis), soc = sum(soc), delta_def = sum(def - soc))]
  
  
}


sim18_ex_fig <- function(){
  
  source("R/libs.R")
  source("R/init.R")
  source("R/util.R")
  
  # CFG
  f_cfgsc <- file.path("./etc/sim18/cfg-sim18-v01.yml")
  l_spec <- config::get(file = f_cfgsc)
  l_spec <- update_sim18_cfg(l_spec)
  l_spec$N_pt <- c(4000, 100, 100)
  l_spec$t_0 <- seq_along(1:sum(l_spec$N_pt)) 
  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$N_pt)
  message("N: ", paste0(cumsum(l_spec$N_pt), collapse = ", "))
  
  # PLOT
  d_cohort <- sim18_cohort(l_spec)$d_cohort
  d_tbl_1 <- copy(d_cohort)
  d_tbl_1[, state := factor(
    state, levels = l_spec$state_opts, labels = names(l_spec$state_opts))]
  d_tbl_1[, trt := factor(trt, levels = l_spec$trt_lab)]
  d_tbl_2 <- d_tbl_1[, .(.N), keyby = .(state, trt, day)]
  d_tbl_2 <- base::merge(
    d_tbl_2, 
    d_tbl_1[, .(N_unit = length(unique(id))), keyby = .(trt)],
    by = "trt", all.x = T
  )
  d_tbl_2[, prop := N/N_unit]
  p_1 <- ggplot(
    d_tbl_2, aes(x = day, y = prop, group = trt, lty = trt)) +
    geom_line(lwd = 0.4) +
    scale_linetype_discrete("") +
    scale_x_continuous("", breaks = seq(1, max(d_tbl_1$day), by = 5)) +
    scale_y_continuous("Proportion", breaks = seq(0, 1, by = 0.1)) +
    facet_wrap(~state, ncol = 1) +
    theme(
      legend.position = "bottom"
    )
  p_2 <- ggplot(d_tbl_1, aes(fill = state, x = day)) +
    geom_bar(position = "fill") +
    scale_fill_discrete("") +
    scale_x_continuous("", breaks = seq(1, max(d_tbl_1$day), by = 5)) +
    scale_y_continuous("", breaks = seq(0, 1, by = 0.1)) +
    facet_wrap(~trt, ncol = 1) +
    theme(
      legend.position = "bottom"
    )
  p_1 + p_2
  
  # TRANSITION PROBS
  t_days <- 1:max(l_spec$visit_days)
  t_gaps <- rep(1, length = max(l_spec$visit_days))
  
  d_tran_trt <- data.table()
  jj <- kk <- 1
  for(jj in seq_along(l_spec$trt_lab)){
    for(kk in seq_along(t_days)){
      
      d_pr <- data.table(
        sim18_transition_matrix(t_days[kk], t_gaps[kk], l_spec$trt_lab[jj], l_spec)
      )
      d_pr[, `:=`(
        from = l_spec$state_lab, 
        trt = l_spec$trt_lab[jj],
        day = t_days[kk])]
      d_pr <- melt(
        d_pr, id.vars = c("trt", "from", "day"), 
        variable.name = "to", value.name = "p_tran"
      )
      d_tran_trt <- rbind(
        d_tran_trt, d_pr
      )
    }
  }
  d_tran_trt[, from := factor(from, levels = l_spec$state_lab)]
  d_tran_trt[, to := factor(to, levels = l_spec$state_lab)]
  d_tran_trt[, trt := factor(trt, levels = l_spec$trt_lab)]
  
  ggplot(
    data = d_tran_trt[day %in% c(1, 4, 7, 14, 21, 28)], 
    aes(x = from, y = to, size = p_tran)) +
    geom_point() +
    scale_x_discrete("Previous State") +
    scale_y_discrete("State") +
    scale_size_continuous("") +
    facet_wrap(trt~day, labeller = label_both, nrow = 3) +
    theme(
      legend.position = "bottom"
    )
  
  
    
  
}


sim18_ex_mod <- function(){
  
  source("R/libs.R")
  source("R/init.R")
  source("R/util.R")
  
  # CFG
  f_cfgsc <- file.path("./etc/sim18/cfg-sim18-v01.yml")
  l_spec <- config::get(file = f_cfgsc)
  l_spec <- update_sim18_cfg(l_spec)
  l_spec$N_pt <- c(4000, 100, 100)
  l_spec$t_0 <- seq_along(1:sum(l_spec$N_pt)) 
  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$N_pt)
  message("N: ", paste0(cumsum(l_spec$N_pt), collapse = ", "))
  
  # MODEL
  l_dat <- sim18_cohort(l_spec)
  d_cohort <- l_dat$d_cohort
  d_obs <- l_dat$d_obs
  m_1 <- cmdstanr::cmdstan_model("stan/sim18-v01.stan")
  l_mod <- sim18_stan_data(d_obs, l_spec)
  # f_1_optim <- m_1$optimize(data = l_mod, jacobian = TRUE)
  # f_1 <- m_1$laplace(data = l_mod, mode = f_1_optim, draws = 2000, refresh = 0)
  f_1 <- m_1$sample(
    l_mod,
    iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
    parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain,
    refresh = 100, show_exceptions = F,
    max_treedepth = 11
  )
  m_post <- f_1$draws(variables = l_spec$non_zero_pars, format = "matrix")
  d_mu_v_tru <- data.table(
    par = l_spec$non_zero_pars,
    tru = c(l_spec$alpha, l_spec$b_trt[-1],  l_spec$b_prev[-1], 
            l_spec$b_time_1, l_spec$b_time_2, 
            # l_spec$b_gap[-1], 
            l_spec$b_prev_time[-1], l_spec$b_trt_time[-1]
    ),
    mu = colMeans(m_post), 
    lo = apply(m_post, 2, function(z){quantile(z, prob = 0.025)}),
    hi = apply(m_post, 2, function(z){quantile(z, prob = 0.975)})
  )
  d_mu_v_tru[, cover := lo < tru & tru < hi]
  kableExtra::kbl(d_mu_v_tru[], format = "simple", digits = 4)
  # 
  
  
  
  
  # print(l_spec)
  # >   print(l_spec)
  # $p_init
  # [1] 0.0 0.4 0.6
  # $desc
  # [1] "Defer at NI Boundary"
  # $nsim
  # [1] 500
  # $mc_cores
  # [1] 40
  # $mcmc_warmup
  # [1] 1000
  # $mcmc_iter
  # [1] 1000
  # $mcmc_chain
  # [1] 1
  # $mcmc_B
  # [1] 1000
  # $nex
  # [1] 10
  # $N_pt
  # [1] 4000  100  100
  # $pt_per_day
  # [1] 1.2
  # $ramp_up_days
  # [1] 60
  # $followup
  # [1] 720
  # $followup_dec
  # [1] 28
  # $visit_days
  # [1]  0  1  2  3  4  5  6  7 14 21 28
  # $state_lab
  # [1] "none"   "mild"   "severe"
  # $trt_active
  # soc  def  dis 
  # TRUE TRUE TRUE 
  # $trt_lab
  # [1] "soc" "def" "dis"
  # $alpha
  # [1] -0.5  1.2
  # $b_trt
  # soc         def         dis 
  # 0.00000000 -0.03682139  0.00000000 
  # $b_prev
  # [1] 0.0 2.2 1.0
  # $b_time_1
  # [1] -0.23
  # $b_time_2
  # [1] 0.003
  # $b_trt_time
  # soc  def  dis 
  # 0.00 0.02 0.00 
  # $b_prev_time
  # [1] 0.00 0.05 0.00
  # $dec_delta_ni
  # [1] 1
  # $dec_thresh_ni
  # [1] 0.985
  # $dec_thresh_fut
  # [1] 0.8
  # $state_opts
  # none   mild severe 
  # 1      2      3 
  # $max_day
  # [1] 28
  # $smry_pars
  # [1] "a"           "b_trt"       "b_prev"      "b_time_1"    "b_time_2"    "b_prev_time" "b_trt_time" 
  # $full_pars
  # [1] "a[1]"           "a[2]"           "b_trt[1]"       "b_trt[2]"       "b_trt[3]"       "b_prev[1]"     
  # [7] "b_prev[2]"      "b_prev[3]"      "b_time_1"       "b_time_2"       "b_prev_time[1]" "b_prev_time[2]"
  # [13] "b_prev_time[3]" "b_trt_time[1]"  "b_trt_time[2]"  "b_trt_time[3]" 
  # $non_zero_pars
  # [1] "a[1]"           "a[2]"           "b_trt[2]"       "b_trt[3]"       "b_prev[2]"      "b_prev[3]"     
  # [7] "b_time_1"       "b_time_2"       "b_prev_time[2]" "b_prev_time[3]" "b_trt_time[2]"  "b_trt_time[3]" 
  # $ex_trial_ix
  # [1]   1  49 105 113 165 179 286 294 464 470
  # $delta_lab
  # [1] "delta_def" "delta_dis"
  # $dur_tru
  # soc      def      dis 
  # 19.79151 18.79142 19.79151 
  # $delta_dur_tru
  # delta_def delta_dis 
  # -1.000094  0.000000 
}


sim18_ex_sim <- function(){
  
  ## SIM SUMMARY
  l_res <- qs::qread(file.path("data/sim18-03/sim18-v01-20260627-175948.qs"))
  l_spec_res <- l_res$l_spec
  d_par <- data.table()
  d_tmp <- copy(l_res$d_post_smry_1[par %in% l_spec_res$non_zero_pars])
  d_tmp[, mu := nafill(mu, type = "locf"), keyby = .(sim, par)]
  
  d_smry <- d_tmp[, .(
    mu = mean(mu),
    lo = quantile(mu, prob = 0.025),
    hi = quantile(mu, prob = 0.975)
  ), keyby = .(ic, par)]
  
  d_smry[, N := cumsum(l_spec_res$N_pt)[ic]]
  d_smry[, stat := sprintf("%.3f (%.3f, %.3f)", mu, lo, hi)]
  d_smry[, par := factor(par, levels = l_spec_res$non_zero_pars)]
  d_tbl <- dcast(d_smry, par ~ N, value.var = "stat")  
  d_tbl[, tru := c(l_spec_res$alpha, l_spec_res$b_trt[-1],  l_spec_res$b_prev[-1], 
                   l_spec_res$b_time_1, l_spec_res$b_time_2, 
                   l_spec_res$b_gap[-1], 
                   l_spec_res$b_prev_time[-1],
                   l_spec_res$b_trt_time[-1]
  )]
  kableExtra::kbl(d_tbl, format = "simple", digits = 4)
}
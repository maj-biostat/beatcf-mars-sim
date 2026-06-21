# uses ordinal model to control the state transitions directly.
# we would need to compute the underlying truth of the target we are trying to 
# estimate.

library(data.table)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(tictoc)

rord_pom <- function(lp, alpha)
{
  p0 <- plogis(alpha[1] - lp)
  p1 <- plogis(alpha[2] - lp) - p0
  p2 <- 1 - p0 - p1
  
  sample(1:3, 1, prob = c(p0, p1, p2)) 
}

sim_cf_trial <- function(
    l_spec
)
{
  
  trt <- sample(l_spec$trt_lab, l_spec$n, replace = TRUE)
  
  dt <- data.table(
    id = rep(seq_len(l_spec$n), each = l_spec$days + 1),
    day = rep(0:l_spec$days, l_spec$n),
    trt = rep(trt, each = l_spec$days + 1)
  )
  
  setorder(dt,id,day)
  
  dt[, state := NA_integer_]
  dt[day==0, state:= sample(l_spec$state_opts, .N, replace = TRUE, prob = l_spec$p_init )]
  
  for(i in seq_len(l_spec$n))
  {
    
    rows <- which(dt$id==i)
    
    for(j in 2:length(rows))
    {
      
      prev <- dt$state[rows[j-1]]
      
      interval <- 1
      
      tt <- dt$day[rows[j]]
      
      trt <- dt$trt[rows[j]]
      
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
        l_spec$b_gap*interval
      
      dt$state[rows[j]] <- rord_pom(lp,l_spec$alpha)
      
    }
    
  }
  
  dt_obs <-
    dt[day %in% l_spec$visit_days]
  
  list(
    daily = dt,
    observed = dt_obs
  )
  
}


transition_matrix <- function(day,
                              trt,
                              l_spec)
{
  
  P <- matrix(0,3,3)
  
  for(prev in l_spec$state_opts){
    
    lp <-
      l_spec$b_trt[trt] +
      l_spec$b_prev[prev] +
      l_spec$b_time_1 * day +
      l_spec$b_time_2 * day^2 +
      l_spec$b_trt_time[trt] * day +
      l_spec$b_gap
    
    p1 <- plogis(l_spec$alpha[1] - lp)
    p2 <- plogis(l_spec$alpha[2] - lp) - p1
    p3 <- 1 - p1 - p2
    
    P[prev,] <- c(p1,p2,p3)
    
  }
  
  # from 
  rownames(P) <- names(l_spec$state_opts)
  # to 
  colnames(P) <- names(l_spec$state_opts)
  P
  
}


state_occupancy <- function(days,
                            l_spec)
{
  
  out <- matrix(NA, days + 1, 3)
  
  out[1,] <- l_spec$p_init
  
  # starting point
  pi <- l_spec$p_init
  
  d_sop <- data.table()
  
  for(trt in l_spec$trt_lab){
    
    for(day in 1:days){
      
      Pt <- transition_matrix(
        day,
        trt,
        l_spec
      )
      
      pi <- drop(pi %*% Pt)
      
      out[day + 1,] <- pi
      
    }
    
    d_sop <- rbind(
      d_sop, 
      data.table(
      trt = trt,
      day = 0:days,
      none = out[,1],
      mild = out[,2],
      severe = out[,3]
    )
    )
  }
  
  
  d_sop
  
  
}


# just set up config
testing_0 <- T
if(testing_0){
  
  l_spec <- list()
  l_spec$n = 600
  l_spec$trt_lab <- c("soc", "def", "dis")
  l_spec$state_opts <- c(none = 1, mild = 2, severe = 3)
  l_spec$days = 28
  l_spec$visit_days = c(0:14,21,28)
  
  l_spec$alpha = c(-0.5,1.2)
  
  l_spec$b_trt = c(
    soc    = 0,
    def    = 0.5,
    dis = 1.2
  )
  
  l_spec$b_prev = c(0, 1.5, 0.1)
  
  # linear term
  l_spec$b_time_1 = -0.4
  # quadratic term
  l_spec$b_time_2 = 0.05
  
  # treatment by time
  l_spec$b_trt_time = c(
    soc = 0,
    def = 0.2,
    dis = 0
  )
  
  # for the data generation 
  l_spec$b_gap = 0
  l_spec$p_init <- c(0, 0.4, 0.6)
}


testing_1 <- F
if(testing_1){
  
  # is unconditioning working?
  
  n_sim <- 10
  
  # truth
  d_sop <- state_occupancy(14, l_spec)
  d_sop <- d_sop[day %in% 1:14]
  d_sop <- melt(d_sop, id.vars = c("trt", "day"), variable.name = "state", value.name = "sop")
  
  
  d_r <- rbindlist(pbapply::pblapply(
    X=1:n_sim, cl = 3, FUN=function(ix) {
      
      sim <- sim_cf_trial(l_spec)
      
      d_tbl_1 <- copy(sim$daily)
      d_tbl_1[, state := factor(state, levels = 1:3, labels = names(l_spec$state_opts))]
      d_tbl_1[, trt := factor(trt, levels = l_spec$trt_lab)]
      
      d_tbl_2 <- d_tbl_1[, .(.N), keyby = .(state, trt, day)]
      d_tbl_2 <- merge(
        d_tbl_2, 
        d_tbl_1[, .(N_unit = length(unique(id))), by = .(trt)],
        by = "trt", all.x = T
      )
      d_tbl_2[, prop := N/N_unit]
      
      d_tbl_2 <- d_tbl_2[day %in% 1:14]
      d_tbl_2
      
    }), idcol = "id_sim")
  
  
  d_smry <- d_r[, .(mu_p = mean(prop)), keyby = .(state, trt, day)]
  d_smry <- merge(
    CJ(state = names(l_spec$state_opts), trt = l_spec$trt_lab, day = 1:14), 
    d_smry, 
    all.x = T, by = c("state", "trt", "day") )
  d_smry[is.na(mu_p), mu_p := 0]
  
  d_smry <- merge(d_smry, d_sop, by = c("trt", "day", "state"), all.x = T)
  
  
  ggplot(d_smry, aes(x = mu_p, y = sop)) + 
    geom_point(size = 0.6) +
    geom_smooth(method = "loess", se = F, lwd = 0.2) +
    facet_grid(trt ~ state)
  #
  
  sim <- sim_cf_trial(l_spec)
  d_fig <- copy(sim$daily)
  d_fig[, state := factor(state, levels = 1:3, labels = names(l_spec$state_opts))]
  d_fig[, trt := factor(trt, levels = l_spec$trt_lab)]
  
  p_1 <- ggplot(d_fig, aes(fill = state, x = day)) +
    geom_bar(position = "fill") +
    scale_fill_discrete("") +
    scale_x_continuous("", breaks = 1:l_spec$days) +
    facet_wrap(~trt, ncol = 1) +
    theme(
      legend.position = "bottom"
    )
  print(p_1)
}


testing_2 <- F
if(testing_2){
  
  # can we recover from model?
  
  sim <- sim_cf_trial(l_spec)
  d_obs <- copy(sim$observed)
  
  setorder(d_obs, id, day)
  
  d_obs[, `:=`(
    prev_state = shift(state, 1L),
    prev_day   = shift(day, 1L)
  ), by = id]
  
  # For the day zero of onset, we have no prev state. We can either assume that they
  # were well the day before or just drop that observation and include it in day 1
  # of follow up. I do the latter.
  d_obs[, gap_len := day - prev_day]
  d_obs <- d_obs[day != 0]
  # gives 16 days of follow up per pt 14 whole days and then 1 obs in wk 3 and 4
  d_obs[, trt_idx := match(trt, l_spec$trt_lab)]
  
  d_obs[, x_time := copy(day)]
  # scaling makes this sample a lot faster but is a pain in the arse
  # for linear predictors with interactions
  # additionally means you need to back scale intercepts
  d_obs[, x_time := scale(x_time)]
  d_obs[, x_gap_time := factor(fifelse(d_obs$gap_len == 1, 1, 2))]
  d_obs[, x_trt := factor(trt_idx)]
  d_obs[, x_prev := factor(prev_state, levels = l_spec$state_opts)]
  
  # b_trt[trt] + b_prev[prev] + 
  #   b_time_1*time + b_time_2*pow(time, 2) +
  #   b_trt_time[trt] .* time +
  #   b_gap[gap_time];
  
  X <- model.matrix(~ x_trt + x_prev +
                      x_time + I(x_time^2) +
                      x_gap_time +
                      x_trt * x_time, 
                    data = d_obs)
  X_mod <- X[, -1]
  
  l_mod <- list(
    N  = nrow(d_obs),
    P = ncol(X_mod),
    X = X_mod,
    y = d_obs$state,
    # back transforms
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
    mu_days = mean(d_obs$day),
    sd_days = sd(d_obs$day),
    T_pred = 28
  )
  l_mod$K = length(l_mod$ix_no_trans)
  
  
  m_1 <- cmdstan_model("./stan/cf-markov-2.stan")
  
  # tic()
  # f_1 <- m_1$sample(
  #   data = l_mod,
  #   chains = 1,
  #   parallel_chains = 1,
  #   iter_warmup = 1000,
  #   iter_sampling = 1000,
  #   refresh = 10,
  #   max_treedepth = 11
  # )
  # toc()
  
  # done in about 20% of the time for the above
  # tic()
  # f_1 <- m_1$pathfinder(
  #   l_mod, num_paths=20,
  #   single_path_draws=200,
  #   history_size=50, max_lbfgs_iters=100,
  #   refresh = 0, draws = 2000
  # )
  # toc()
  
  tic()
  f_1_optim <- m_1$optimize(data = l_mod, jacobian = TRUE)
  f_1 <- m_1$laplace(data = l_mod, mode = f_1_optim, draws = 2000, refresh = 0)
  toc()
  
  smry_pars <- c("a", "b_trt", "b_prev", "b_time_1", "b_time_2", "b_gap", "b_trt_time")
  d_post <- data.table(
    f_1$draws(variables = smry_pars, format = "matrix")
  )
  d_smry <- melt(d_post, measure.vars = names(d_post))
  d_smry <- d_smry[, .(
    mu = mean(value),
    q_025 = quantile(value, prob = 0.025),
    q_975 = quantile(value, prob = 0.975)
  ), keyby = variable]
  
  par_tru <- c(
    l_spec$alpha, 
    l_spec$b_trt,
    l_spec$b_prev,
    l_spec$b_time_1,
    l_spec$b_time_2,
    l_spec$b_gap,
    l_spec$b_trt_time
  )
  names(par_tru) <- d_smry$variable
  
  d_smry[, truth := par_tru[variable]]
  d_smry[, cover := (q_025 < truth & truth < q_975) | mu == truth]
  d_smry[]
    
  
  # generate transition matrix for every draw
  
  # the gap parameter would probably be excluded for daily
  
  # obtain the sop
    
  
  
  
  
  
}



library(data.table)






get_sim02_trial_data <- function(
    l_spec
){
  
  # cohort
  if(is.null(l_spec$ic)){
    ic <- 1
  } else {
    ic <- l_spec$ic
  }
  if(is.null(l_spec$is)){
    is <- 1
  } else {
    is <- l_spec$is
  }
  if(is.null(l_spec$ie)){
    ie <- 1
  } else {
    ie <- l_spec$ie
  }
  if(is.null(l_spec$t0)){
    t0 <- 1
  } else {
    t0 <- l_spec$t0
  }
  if(is.null(l_spec$tfu)){
    tfu <- 1
  } else {
    tfu <- l_spec$tfu
  }
  
  trt_opts <- 1:length(l_spec$p_trt_alloc)
  
  sigmas <- c(l_spec$b_s0, l_spec$b_s1)
  s <- diag(sigmas)
  r <- matrix(c(1, l_spec$b_rho, l_spec$b_rho, 1), nrow = 2)
  sigma <- s %*% r %*% s
  
  d_ref <- data.table(
    # interim id
    ic = ic, 
    # unit id
    id = is:ie,
    trt = sample(trt_opts, size = l_spec$N[ic], replace = T, prob = l_spec$p_trt_alloc),
    t0 = t0,
    tfu = tfu,
    MASS::mvrnorm(l_spec$N[ic], mu = c(0, 0), Sigma = sigma))
  
  setnames(d_ref, c("V1", "V2"), c("u_0", "u_1"))
  
  # unit specific int and slope over 1 year.
  # but i don't think you could meaningfully extrapolate a 
  # year on year impact?
  d_ref[, b_0 := l_spec$b_0 + u_0]
  d_ref[, b_1 := l_spec$b_time + u_1 + l_spec$b_trt[trt]]
  
  fu <- c(0, 1)
  d <- d_ref[, .(fu = fu), by = .(ic, id, trt, t0, tfu, b_0, b_1)]
  
  d[, ia := NA_integer_]
  d[, t_anlys := NA_real_]
  
  d[, y := b_0 + b_1 * fu + rnorm(.N, 0, l_spec$b_se)]
  
  d[, y_mis := rbinom(.N, 1, l_spec$pr_ymis)]
  
  setcolorder(d, c("ic", "ia", "id", "trt", "t0", "tfu", "t_anlys"))
  
  d
}



get_sim02_stan_data <- function(d_mod){
  
  # all pts with any missingness are dropped
  dd <- dcast(
    d_mod[!(id %in% d_mod[y_mis == 1, id])], 
    ic + id + t0 + trt ~ fu, value.var = "y")
  setnames(dd, c("0", "1"), c("y_pre", "y"))
  
  # mean centre pre so that the intercept makes a bit more sense.
  dd[, y_pre := y_pre - mean(y_pre)]
  
  ld <- list(
    N = nrow(dd),
    y = dd$y,
    y_pre = dd$y_pre,
    trt = dd$trt,
    prior_only = 0
  )
  
  list(
    dd = dd,
    ld = ld
  )
  
}

prototype_2_data <- function(){
  
  library(data.table)
  library(ggplot2)
  library(ggh4x)
  library(MASS)
  library(sn)
  library(pbapply)
  library(gt)
  library(gtsummary)
  library(pander)
  library(cmdstanr)
  
  
  
  get_dummy_dat <- function(){
    
    
    # population baseline and trend
    a_0     <- 80      # 
    a_1     <- -1
    u_0     <- 0
    u_1     <- 0
    s_0 <- 0.25      
    s_1 <- 0.05      
    rho    <- 0.2   
    
    # both trt strategies make fev a little bit worse
    b_trt     <- c(0, -0.4, -0.4)
    
    # std dev within participants
    s_e <- 0.6  
    
    # data
    mu     <- c(u_0, u_1)          # combine the means in a vector
    sigmas <- c(s_0, s_1)  # combine the std devs in a vector
    
    s <- diag(sigmas)      # standard deviation matrix
    # correlation matrix
    r <- matrix(c(1, rho,  rho, 1), nrow = 2)
    
    # now matrix multiply s and r to get a covariance matrix
    sigma <- s %*% r %*% s
    
    # how many participants would you like?
    N <- 100
    
    # unit specific baseline and slope that characterises the 
    # trajectory of fev over time, irrespective of trt
    # i.e. people are expected to decline over time.
    d_ref <- data.table(
      id = 1:N, 
      trt = sample(1:3, size = N, replace = T),
      MASS::mvrnorm(N, mu = mu, Sigma = sigma))
    setnames(d_ref, c("V1", "V2"), c("u_0", "u_1"))
    
    # unit specific int and slope over 1 year.
    # but i don't think you could meaningfully extrapolate a 
    # year on year impact?
    d_ref[, b_0 := a_0 + u_0]
    d_ref[, b_1 := a_1 + u_1 + b_trt[trt]]
    
    times <- c(0, 1)
    d <- d_ref[, .(t = times), by = .(id, trt, b_0, b_1)]
    d[, y := b_0 + b_1 * t + rnorm(.N, 0, s_e)]
    
    d
  }
  
  
  
  
  d <- get_dummy_dat()

  d[, .(mu = mean(y)), keyby = .(trt, t)]

  min_y <- min(floor(d$y))
  max_y <- max(ceiling(d$y))

  d_fig <- copy(d)
  d_fig[, trt := factor(trt)]

  ggplot(d, aes(x = t, y = y, group = id)) +
    geom_line() +
    # geom_smooth(
    #   data = d, aes(x = t, y = y), inherit.aes = F,
    #   se = F) +
    scale_y_continuous("", breaks = seq(min_y, max_y, by = 1)) +
    ggh4x::facet_wrap2(~trt)
  #
  
  
  # assume treatment is randomized at baseline to individuals 
  # sampled from the same population
  # run a standard pre/post (ancova) analysis
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  m1 <- cmdstanr::cmdstan_model("stan/sim02-v01.stan")
  i <- 1
  r <- pblapply(1:100, FUN = function(i){

    d <- get_dummy_dat()
    d_tmp <- dcast(d, id + trt ~ t, value.var = "y")
    setnames(d_tmp, c("0", "1"), c("y_pre", "y"))
    # d_tmp[, delta := y - y_pre]
    # d_tmp[, mean(delta), keyby = trt]
    
    ld <- list(
      N = nrow(d_tmp),
      y = d_tmp$y,
      y_pre = d_tmp$y_pre,
      trt = d_tmp$trt,
      pri_a = c(80, 10),
      pri_b_trt = c(0, 10),
      pri_b_pre = c(0, 10),
      prior_only = 0
    )
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", i)
    
    # fit model - does it matter that I continue to fit the model after the
    # decision is made...?
    
    f_1 <- m1$sample(
      ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = output_dir_mcmc,
      output_basename = foutname
    )
    
    # f_1 <- m1$pathfinder(
    #   ld, num_paths=20, 
    #   single_path_draws=200,
    #   history_size=50, max_lbfgs_iters=100,
    #   refresh = 0, draws = 2000)
    
    d_post <- data.table(f_1$draws(
      variables = c("mu"), format = "matrix"))
    d_fig <- melt(d_post, measure.vars = names(d_post))
    ggplot(d_fig, aes(x = value, group = variable, col = variable)) +
      geom_density()
    
    d_post <- data.table(f_1$draws(
      variables = c("delta_2_1", "delta_3_1"), format = "matrix"))
    d_fig <- melt(d_post, measure.vars = names(d_post))
    ggplot(d_fig, aes(x = value, group = variable, col = variable)) +
      geom_density()
    
    # diffs in means between trt and soc
    # we expect mu2 and mu3 to be a little worse than mu1 
    # but not materially so
    
    # probability that trt2 is no worse than ni margin relative to soc
    d_post[, .(pr_2_1 = mean(delta_2_1 > -0.5),
               pr_3_1 = mean(delta_3_1 > -0.5)) ]
    
  }, cl = 4)
  
  d_res <- rbindlist(r)
  d_res[, `:=`(win_2 = pr_2_1 > 0.95, win_3 = pr_3_1 > 0.95)]
  colMeans(d_res[, .(win_2, win_3)])
#  
    
}

prototype_1_data <- function(){
  
  library(data.table)
  library(ggplot2)
  library(ggh4x)
  library(MASS)
  library(sn)
  library(pbapply)
  library(gt)
  library(gtsummary)
  library(pander)
  
  # y_{i} \sim normal(\mu_{i}, \sigma_e)
  # \mu_{i} = a_0 + 
    # b_1 * I(18<age<30) + b_2 * I(age>=30) +
    # b_3 * I(lung<70) + 
    # b_4 * I(pseudo) +
    # b_5 * I(pseudo & lung<70) 
    # ...?
    
  # assume:
  # 30% of pop is <18. 
  # of this cohort 3% have lung func <70% and 10% have pseudo
  
  # 40% of pop is 18<=age<30
  # of this cohort 15% have lung func <70% and 30% have pseudo
  
  # 30% of pop is >=30
  # of this cohort 20% have lung func <70% and 40% have pseudo
  
  # what is the population size of the cohort we are interested in?
  N <- 5000
  d <- data.table(
    id = 1:N  
  )
  
  # want some kind of skewed distribution for age, solve based on gamma
  # Target quantiles
  q_probs <- c(0.3, 0.7)       # 30% less than 18, 70% less than 30
  tau <- c(18, 30)        # age cut-points
  
  # compute distance from target quantiles for given gamma params
  # objective <- function(par) {
  #   shape <- par[1]
  #   scale <- par[2]
  #   q_target <- qgamma(q_probs, shape = shape, scale = scale)
  #   sum((q_target - tau)^2)  # sum of squared errors
  # }
  # 
  # # Find gamma shape & scale that match targets
  # fit <- optim(par = c(5, 2), fn = objective, method = "L-BFGS-B",
  #              lower = c(0.1, 0.1))
  
  shape_est <- 4.501337 # fit$par[1]
  scale_est <- 5.628841 # fit$par[2]
  
  # simulate age based on pars (truncated so that silly ages are avoide)
  age_max <- 65
  d[, age_0 := rgamma(.N, shape = shape_est, scale = scale_est)]
  # truncate
  d <- d[age_0 < age_max]
 
  d[age_0 < 18, pseudo_0 := rbinom(.N, 1, 0.20)]
  d[age_0 >= 18 & age_0 < 30, pseudo_0 := rbinom(.N, 1, 0.30)]  
  d[age_0 >= 30, pseudo_0 := rbinom(.N, 1, 0.45)]  
  
  d[, age_cat := fcase(
    age_0 < 18, 1L,
    age_0 >= 18 & age_0 < 30, 2L,
    age_0 >= 30, 3L
  )]
  
  # reference group is:
  # below 18 (age 0)
  # lung func > 70
  # no pseudo
  # say ppfev1 around 97%
  
  b_0 <- 97
  
  # initial annual decline
  b_age <- -0.75
  # delta_age[1] change in annual decline for 18 <= age < 30 relative to b_age_0
  # delta_age[2] change in annual decline for 18 <= age < 30 relative to b_age_0 + b_age_1
  delta_age_1 <- 0.35
  delta_age_2 <- 0.32
  
  # assume the effect of pseudo is invariant to age
  b_pseudo <- -4
  
  # pmax(0, age - tau[1]) etc are the basis funcs for segmented regn
  
  d[, mu_0 := b_0 + 
      b_age * age_0 + 
      delta_age_1 * pmax(0, age_0 - tau[1]) + 
      delta_age_2 * pmax(0, age_0 - tau[2]) +
      b_pseudo * pseudo_0
      ]
  
  # assume homogeneous variance even though this probably isn't very likely
  # but ...
  d[, y_0 := rnorm(.N, mu_0, 10)]
  
  
  # followup to 6 mnths
  d[, age_6 := age_0 + 0.5]
  d[, pseudo_6 := fifelse(
    age_6 < 18 & pseudo_0 == 1,
    rbinom(.N, 1, 0.4),
    rbinom(.N, 1, 0.2)
  )]
  d[, pseudo_6 := fifelse(
    age_6 >= 18 & age_6 < 30 & pseudo_0 == 1,
    rbinom(.N, 1, 0.45),
    rbinom(.N, 1, 0.3)
  )]
  d[, pseudo_6 := fifelse(
    age_6 >= 18 & age_6 < 30 & pseudo_0 == 1,
    rbinom(.N, 1, 0.45),
    rbinom(.N, 1, 0.45)
  )]
  
  d[, mu_6 := mu_0 + 
      b_age * age_6 + 
      delta_age_1 * pmax(0, age_6 - tau[1]) + 
      delta_age_2 * pmax(0, age_6 - tau[2]) +
      b_pseudo * pseudo_6
  ]
  
  d[, y_6 := rnorm(.N, mu_6, 10)]
  
  # followup to 6 mnths
  d[, age_12 := age_0 + 1]
  d[, pseudo_12 := fifelse(
    age_12 < 18 & pseudo_6 == 1,
    rbinom(.N, 1, 0.4),
    rbinom(.N, 1, 0.2)
  )]
  d[, pseudo_12 := fifelse(
    age_12 >= 18 & age_12 < 30 & pseudo_6 == 1,
    rbinom(.N, 1, 0.45),
    rbinom(.N, 1, 0.3)
  )]
  d[, pseudo_6 := fifelse(
    age_12 >= 18 & age_12 < 30 & pseudo_6 == 1,
    rbinom(.N, 1, 0.45),
    rbinom(.N, 1, 0.45)
  )]
  
  d[, mu_12 := mu_0 + 
      b_age * age_12 + 
      delta_age_1 * pmax(0, age_12 - tau[1]) + 
      delta_age_2 * pmax(0, age_12 - tau[2]) +
      b_pseudo * pseudo_12
  ]
  
  d[, y_12 := rnorm(.N, mu_12, 10)]
  
  
  d_fig <- copy(
    d[, .(id, 
          age_0, age_6, age_12, 
          pseudo_0, pseudo_6, pseudo_12,
          y_0, y_6, y_12
          )])
  d_fig[, `:=`(pseudo_0 = as.numeric(pseudo_0), 
               pseudo_6 = as.numeric(pseudo_6), 
               pseudo_12 = as.numeric(pseudo_12))]
  d_fig <- melt(
    d_fig, 
    id.vars = c("id"),
    measure.vars = patterns(
      age = "^age_",
      pseudo = "^pseudo_",
      y = "^y_",
      cols=names(d_fig)
    ),
    variable.name = "t")
  
  setkey(d_fig, id)
  
  d_fig[, age_cat := fcase(
    age < 18, 1L,
    age >= 18 & age < 30, 2L,
    age >= 30, 3L
  )]
  
  ggplot(d_fig[id %in% 1:100], aes(x = t, y = y, group = id)) +
    geom_line(size = 0.3)  +
    facet_wrap2(vars(age_cat), labeller = label_both)
  
  
  d_tbl <- copy(d)
  d_tbl[, age_cat := factor(age_cat, labels = c(
    "[0,18)", "[18,30)","[30,65]"
  ))]
  d_tbl[, pseudo := factor(pseudo, labels = c("no-pseudo", "pseudo"))]
  setnames(d_tbl, "y", "ppFev1")
  d_tbl |>
    tbl_summary(
      include = c(age_cat, ppFev1),
      by = pseudo,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})"
      )
      ) |>
    add_n() |>
    add_overall() |>
    modify_spanning_header(c("stat_1", "stat_2") ~ "**Pseudomonas**") |>
    modify_header(label = "Variable")
  #
  d_tbl |>
    tbl_summary(
      include = c(pseudo, ppFev1),
      by = age_cat,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})"
      )
    ) |>
    add_n() |>
    add_overall() |>
    modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Age**") |>
    modify_header(label = "Variable")
  #
  d_tbl |>
    tbl_strata(
      strata = age_cat,
      .tbl_fun =
        ~ .x |>
        tbl_summary(
          include = c(ppFev1),
          by = pseudo, missing = "no") |>
        add_n() ,
      .header = "**{strata}**, N = {n}"
    )  
  
  
  
  
  # alternatively
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # So now consider unit level progression based on their present
  # (baseline) state. Perhaps consider multiple timepoints over a 
  # 1 year with the rate of decline depending on age group
  
  # y_{ti} \sim normal(\mu_{ti}, \sigma_e)
  # \mu_{ti} = \mu_0 + 
  #     g_time_age_1 * time_1 * I(age_cat == 1) +
  #     g_time_age_1 * time_1 * I(age_cat == 2) +
  #     g_time_age_1 * time_1 * I(age_cat == 3) +
  #     
  #     delta_time_age_1 * time_1 +
  #     g_time_age_1 * time_1 +
  #     delta_age_1 * pmax(0, age - tau[1]) + 
    delta_age_2 * pmax(0, age - tau[2]) +
  #     u_{0i} + u_{1i} time_{ti} +
  #     
  # [u_{0i}] \sim MVN([0], \Sigma)
  # [u_{1i}]          [0]
  # \Sigma = SRS
  # S = [\sigma_0, 0]
  #     [0,        \sigma_1]
  # R = [1,   \rho]
  #     [\rho, 1]
  
  # where \beta_0 is the population-level intercept (initial status) and 
  # \beta_1  is the population-level slope (change over time). The 
  # u_{0i} and u_{1i} terms are the participant-level deviations around the 
  # population-level intercept and slope. 
  
  # \sigma_0 and \sigma_1 capture differences between participants in 
  # their intercepts and slopes, whereas \sigma_e captures the differences 
  # within participants over time that occur apart from their linear trajectories.
  
  
  b0     <- 0      # starting point (average intercept)
  b1     <- 1      # growth over time (average slope)
  sigma0 <- 1      # std dev in intercepts
  sigma1 <- 1      # std dev in slopes
  rho    <- -0.5    # correlation between intercepts and slopes
  sigma_e <- 0.75  # std dev within participants
  
  
  mu     <- c(b0, b1)          # combine the means in a vector
  sigmas <- c(sigma0, sigma1)  # combine the std devs in a vector
  
  s <- diag(sigmas)      # standard deviation matrix
  r <- matrix(c(1, rho,  # correlation matrix
                rho, 1), nrow = 2)
  
  # now matrix multiply s and r to get a covariance matrix
  sigma <- s %*% r %*% s
  
  # how many participants would you like?
  N <- 100
  
  # make the simulation reproducible
  set.seed(1)
  
  # unit level random intercept and slope
  d_ref <- data.table(
    id = 1:N,
    MASS::mvrnorm(n_id, mu = mu, Sigma = sigma)  
  )
  setnames(d_ref, c("V1", "V2"), c("u_0", "u_1"))
  
  # timepoints
  t <- 3
  
  # data set
  d <- d_ref[rep(1:.N, each = t)]
  d[, y := rnorm(.N, u_0 + u_1)]
  
  
}

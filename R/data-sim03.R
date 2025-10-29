
if(exists("prefix_r")){
  source(paste0(prefix_r, "/libs.R"))
  source(paste0(prefix_r, "/data.R"))
} else {
  source("./R/libs.R")
  source("./R/data.R")
}


get_sim03_trial_data <- function(
    l_spec
){
  # cohort
  if(is.null(l_spec$ic)){ ic <- 1 } else { ic <- l_spec$ic }
  if(is.null(l_spec$is)){ is <- 1 } else { is <- l_spec$is }
  if(is.null(l_spec$ie)){ ie <- 1 } else { ie <- l_spec$ie }
  if(is.null(l_spec$t0)){ t0 <- 1 } else { t0 <- l_spec$t0 }
  
  trt_opts <- 1:length(l_spec$p_trt_alloc)
  
  N_ic <- length(is:ie)
  
  d_i <- data.table(
    # interim id
    ic = ic, 
    # unit id
    id = 1:N_ic,
    t0 = t0
  )
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    l_spec$age_lwr, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig, lower.tail = T)
  p_gt_age_upr <- plnorm(
    l_spec$age_upr, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(N_ic, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  baseline_age <- qlnorm(u, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig)
  d_i[, age0 := baseline_age]
  trt_alloc <- trt_opts[rep(trt_opts[as.logical(l_spec$p_trt_alloc)], len = N_ic)]
  d_i[, trt := trt_alloc]
  d_i[, gamma_i := rnorm(.N, 0, l_spec$sigma_i)]
  # latent baseline measure on which the data is generated
  d_i[, mu0 := l_spec$b_0 + l_spec$b_age * log(age0) + gamma_i]
  # what we would actually observe and therefore model
  d_i[, y0 := rnorm(.N, mu0, l_spec$sigma)]
  
  d <- CJ(
    # unit id
    id = 1:N_ic,
    # spirometry observation point - excludes baseline
    t_id = seq_along(l_spec$t_sprty_obs) 
  )
  d[, t_obs := l_spec$t_sprty_obs[t_id]]
  d[, t_fu := t0 + 365*t_obs]
  d <- d[d_i, on = "id"]
  
  
  gamma_trt <- l_spec$b_trt[cbind(d$trt, d$t_id)]
  d[, b_trt := gamma_trt]
  d[, b_time := l_spec$b_time[t_id]]
  d[, mu := mu0 + b_time + b_trt]
  
  
  # Covariance matrix for AR(1) structure only applies to the follow up
  Sigma <- outer(
    1:(length(l_spec$t_sprty_obs)), 
    1:(length(l_spec$t_sprty_obs)), 
    function(i, j){
      l_spec$sigma^2 * l_spec$rho^abs(i - j)
  } )
  
  # outcome
  # i <- 600
  for(i in 1:N_ic){
    
    d[id == i, y := as.numeric(rmvnorm(n = 1, mean = d[id == i, mu], sigma = Sigma))]
  }
  
  # missing
  d[, y_mis := rbinom(.N, 1, l_spec$pr_ymis)]
  
  d[, ia := NA_integer_]
  d[, t_anlys := NA_real_]
  d[, id := rep(is:ie, each = length(l_spec$t_sprty_obs))]
  
  setcolorder(
    d,
    c("ic", "id", "t_id", "t_obs", "t0", "t_fu", "age0", 
      "trt", "gamma_i", "mu0", "y0"))
  
  d
}





get_sim03_stan_data <- function(d_mod){
  
  # mean centre pre so that the intercept makes a bit more sense.
  dd <- copy(d_mod)
  dd[, trt := factor(trt)]
  
  X <- model.matrix(~y_pre + trt, data = dd)
  
  ld <- list(
    N = nrow(dd),
    y = dd$y_post,
    Kx = ncol(X) - 1,
    X = X[, 2:ncol(X)],
    prior_only = 0
  )
  
  list(
    dd = dd,
    ld = ld
  )
  
}

get_prototype_cfg <- function(){
  
  l_spec <- list()
  
  l_spec$desc <- "Test"
  l_spec$nsim <- 1
  l_spec$nex <- 1
  
  # N by analysis
  l_spec$N <- 600
  
  l_spec$pt_per_day <-  0.57
  l_spec$ramp_up_days <-  120
  
  l_spec$fu_days <- 365
  l_spec$fu_days_lwr <- 351
  l_spec$fu_days_upr <- 379
  
  # missingness
  l_spec$pr_ymis <- 0.1
  
  # trt alloc
  l_spec$p_trt_alloc <- rep(1, 3)
  
  # model parameters ppfev
  # reference level at age 0
  l_spec$b_0 <- 120
  # backgroun temporal change reduces the population level ppfev1
  l_spec$b_time <- c(-0.25, -0.5, -0.75, -1)
  # coefficient on log-age
  l_spec$b_age <- -14
  # l_spec$b_ln_age <- -2
  # treatment effects
  l_spec$b_trt <- rbind(
    c(0, 0, 0, 0),
    c(1, 2, 3, 4),
    c(-1, -2, -3, -4)
  )
  
  
  # age class distribution
  # quantile_age_data <- function(age = 20, mu, sig, lower.tail = T){
  #   p_gt_age <- plnorm(age, mu, sig, lower.tail)
  #   p_gt_age
  # } 
  # plot_age_data <- function(mu = 1, sig = 1){
  #   x <- seq(0, 100, len = 1000)
  #   y <- dlnorm(x, mu, sig)
  #   plot(x, y, type = "l", xlim=c(0, 100), xaxt='n' )
  #   abline(v = 5, col = "red", lwd = 0.3)
  #   abline(v = 10, col = "red", lwd = 0.3)
  #   abline(v = 15, col = "red", lwd = 0.3)
  #   axis(side = 1, at = c(seq(0, 20, by = 1), seq(25, 100, by = 5)) )
  # }
  # plot_age_data(mu = 2.5, sig = 0.75)
  # quantile_age_data(c(10, 50, 90), 2.5, 0.75, lower.tail = F)

  # baseline log age/ log sig (log normally distributed age)
  # l_spec$age_mu <- -1.8
  # l_spec$age_sig <- 0.6
  # l_spec$age_lwr <- 0.03
  # l_spec$age_upr <- 0.85
  
  l_spec$age_mu <- 2.5
  l_spec$age_sig <- 0.75
  l_spec$age_lwr <- 3
  l_spec$age_upr <- 85
  # l_spec$age_ref <- 25
   
  # baseline + quarterly spirometry (0, 3, 6, 9, 12 months)
  l_spec$t_sprty_obs <- c(0.25, 0.5, 0.75, 1)
  
  l_spec$sigma_i <- 2
  
  # AR(1) correlation structure for errors
  l_spec$rho <- 0.85      # autocorrelation between adjacent time points
  l_spec$sigma <- 6       # standard deviation at each time point
  
  # priors on log-odds/log-or
  l_spec$prior <- list()
  
  l_spec$prior$pri_b_0 <- c(100, 7)
  l_spec$prior$pri_b_mu <- rep(0, 3)
  l_spec$prior$pri_b_sd <- rep(10, 3)

  
  # decision hurdle
  l_spec$delta <- list()
  l_spec$delta$ni <- -0.5
  l_spec$delta$inf <- -0.5
  
  # evidentiary requirement
  l_spec$thresh <- list()
  l_spec$thresh$ni <- 0.98
  l_spec$thresh$inf <- 0.7
  
  lambda = l_spec$pt_per_day
  # ramp up over x months 
  rho = function(t) pmin(t/l_spec$ramp_up_days, 1)
  
  # lambda = 0.57
  # # ramp up over x months 
  # rho = function(t) pmin(t/120, 1)
  # 
  # rr <- unlist(pblapply(1:100, cl = 4, FUN=function(ii){
  #   ttt <- get_enrol_time(sum(l_spec$N), lambda, rho)
  #   max(ttt)
  # }))
  # mean(rr) / 365
  
  # day of enrolment
  loc_t0 <- get_enrol_time(sum(l_spec$N), lambda, rho)
  # follow up time
  loc_t <- loc_t0 + runif(
    length(loc_t0), l_spec$fu_days_lwr, l_spec$fu_days_upr)
  
  l_spec$ic <- 1 # interim number
  
  # next chunk of data on pts.
  if(l_spec$ic == 1){
    # starting pt index in data
    l_spec$is <- 1
    l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
  } else {
    l_spec$is <- l_spec$ie + 1
    l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
  }
  
  # study date and follow up date
  l_spec$t0 <- loc_t0[l_spec$is:l_spec$ie]
  l_spec$tfu <- loc_t[l_spec$is:l_spec$ie]
  
  l_spec
}



prototype_1_data <- function(){
  
  l_spec <- get_prototype_cfg()
  str(l_spec)

  d <- get_sim03_trial_data(l_spec)
  d[, trt := factor(trt)]
  d[, t_obs := factor(t_obs)]
  
  # f_1_f <- nlme::gls(y ~ log(age0) + trt + t_obs, d_mod2,
  #                    correlation = nlme::corAR1(form = ~ 1 | id))
  
  d[, y0 := scale(y0)]
  f_1_f <- nlme::gls(y ~ y0 + log(age0) + t_obs*trt, d,
                     correlation = nlme::corAR1(form = ~ 1 | id))
  
  s <- summary(f_1_f)

  # plot_1(d)
  
  X <- model.matrix(~ y0 + log(age0) + t_obs*trt, data = d)
  
  # assume treatment is randomized at baseline to individuals 
  # sampled from the same population
  # run a standard pre/post (ancova) analysis
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  
  m1 <- cmdstanr::cmdstan_model("stan/sim03-v06b.stan")  
  
  # // baseline (y0, i.e. first observation of ppFEV1)
  # // log(age0)     
  # // t_obs0.25      
  # // t_obs0.5       
  # // t_obs0.75      
  # // t_obs1         
  # // trt2           
  # // trt3           
  # // t_obs0.25:trt2
  # // t_obs0.5:trt2  
  # // t_obs0.75:trt2 
  # // t_obs1:trt2    
  # // t_obs0.25:trt3
  # // t_obs0.5:trt3 
  # // t_obs0.75:trt3 
  # // t_obs1:trt3   
  # 
  # X <- model.matrix(~ y0 + log(age0) + )
  
  # goal to simplify the data generation and then get the stan model to
  # recover parameters - compare to gls
  
  # all those that have completed followup to 1 year
  lsd <- list(
    N = d[, .N],
    J = length(l_spec$t_sprty_obs)-1,
    P = ncol(X)-1,
    X = X[, -1],
    y = d[, y],
    N_1 = d[t_id == 1, .N],
    N_2 = d[t_id == 2, .N],
    N_3 = d[t_id == 3, .N],
    N_4 = d[t_id == 4, .N],
    ix_1 = d[t_id == 1, which = T],
    ix_2 = d[t_id == 2, which = T],
    ix_3 = d[t_id == 3, which = T],
    ix_4 = d[t_id == 4, which = T]
  )
  
  foutname_1 <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"),
    "-sim-", 1)
  
  f_1 <- m1$sample(
    lsd, iter_warmup = 1000, iter_sampling = 1000,
    parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = T,
    max_treedepth = 10,
    output_dir = output_dir_mcmc,
    output_basename = foutname_1,
    init = list(
      list(
        rho = 0.5,
        sigma = 1)
    ) 
  )
  cbind(s$coefficients, f_1$summary(variables = c("b0", "b")))
  f_1$summary(variables = c("sigma", "rho"))
  s$sigma
  
  
  # f_1_optim <- m1$optimize(data = ld, jacobian = TRUE)
  # f_1_a <- m1$laplace(data = ld, mode = f_1_optim, draws = 2000)
  # f_1_a$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho", "nu"))

  # f_1_b <- m1$variational(data = ld)
  # f_1_b$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho", "nu"))
  # 
  # f_1_c <- m1$pathfinder(data = ld)
  # f_1_c$output()
  # f_1_c$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho"))
  
  
  
  d_post_1 <- data.table(f_1$draws(variables = c("b0", "b"),  format = "matrix"))
  
  d_fig_1 <- melt(d_post_1, measure.vars = names(d_post_1))
  
  ggplot(d_fig_1, aes(x = value)) +
    geom_density() +
    ggh4x::facet_wrap2(~variable, ncol = 2, scales = "free_x")
}


plot_1 <- function(d){
  
  # d[, .(mu = mean(y)), keyby = .(trt, t_obs)]
  # min_y <- min(floor(d$y))
  # max_y <- max(ceiling(d$y))
  
  d_fig <- copy(d)
  d_fig[, trt := as.numeric(trt)]
  d_fig[, t_obs := as.numeric(t_obs)]
  
  d_fig[, age_cat := cut(age0, breaks  = c(0, 0.1, 0.2, 0.3, 0.5, 0.8))]
  d_fig[, age := age0 + t_obs]
  
  rnd_id <- sample(unique(d$id), size = 50)
  p1 <- ggplot(d_fig[id %in% rnd_id], aes(x = t_obs, y = y)) +
    geom_line(aes(group = id), lwd = 0.2) +
    geom_point(size = 0.5) +
    # scale_x_continuous(lim = c(0, 90), breaks = seq(0, 90, by = 5)) +
    # scale_y_continuous(lim = c(25, 130), breaks = seq(25, 125, by = 25)) +
    ggh4x::facet_wrap2(~trt, ncol = 3)
  #
  
  p2 <- ggplot(d_fig, aes(x = age, y = y)) +
    geom_line(aes(group = id), lwd = 0.2) +
    # scale_x_continuous(lim = c(0, 90), breaks = seq(0, 90, by = 5)) +
    # scale_y_continuous(lim = c(25, 130), breaks = seq(25, 125, by = 25)) +
    ggh4x::facet_wrap2(~trt, ncol = 1)
  #
  p2 + p1
  #
  
  
}

prototype_2_data <- function(){
  
  # expt with dealing with the data in long format
  # so can work towards handling missing data.
  
  l_spec <- get_prototype_cfg()
  str(l_spec)
  
  d <- get_sim03_trial_data(l_spec)
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  
  m1 <- cmdstanr::cmdstan_model("stan/sim03-v07.stan")  
  
  d_mod <- copy(d)
  d_mod[, trt := factor(trt, levels = 1:3)]
  X <- model.matrix(~ log(age0) + t_obs + trt, data = d_mod)[, -1]
  
  d_ix <- d_mod[, .(ix_s = .I[1], ix_e = .I[.N]), by=id]
  
  ld <- list(
    N = d_mod[, .N],
    Q = length(unique(d_mod$id)),
    y = d_mod$y,
    s_id = d_ix$ix_s,
    e_id = d_ix$ix_e,
    X = X,
    P = ncol(X)
  )
  foutname_1 <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"),
    "-sim-", 1)
  
  f_1 <- m1$sample(
    ld, iter_warmup = 1000, iter_sampling = 1000,
    parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = T,
    max_treedepth = 10,
    output_dir = output_dir_mcmc,
    output_basename = foutname_1
  )
  f_1$summary(
    variables = c("b0", "b", "sigma", "rho"))
  
  
  #
  
  
}

tmp <- function(){
  library(data.table)
  set.seed(126)
  n = 10
  a <- 20  #intercept
  b <- 0.2  #slope
  x <- round(runif(n, 1, n), 1)  #values of the year covariate
  year <- 1:n
  sigma <- 20
  rho <- 0.8
  
  library(nlme)
  ## define a constructor for a first-order
  ## correlation structure
  ar1 <- corAR1(form = ~year, value = rho)
  str(ar1)
  ## initialize this constructor against our data
  AR1 <- Initialize(ar1, data = data.frame(year))
  str(AR1)
  ## generate a correlation matrix
  V <- corMatrix(AR1)
  V
  
  ## Cholesky factorization of V
  Cv <- chol(V)
  ## simulate AR1 errors
  e <- t(Cv) %*% rnorm(n, 0, sigma)  # cov(e) = V * sig^2
  ## generate response
  y <- a + b * x + e
  data.temporalCor = data.table(y = y, x = x, year = year)
              
    
  
  
  rho <- 0.7
  d_1 <- data.table()
  sig_gamma <- 2
  # sd marginal
  sig_marg <- sqrt((sig_gamma^2) / (1 - rho^2))
  
  for(i in 1:1e4){
    
    n_ptcl <- 100
    g <- numeric(n_ptcl)
    
    g[1] <- rnorm(1, 0, sig_marg)
    
    for(j in 2:n_ptcl){
      g[j] = rho*g[j-1] + rnorm(1, 0, sig_gamma)
    }
    d_1 <- rbind(d_1, data.table(seq = i, g = g))
  }
  d_1[, var(g), keyby = seq][, mean(V1)]
  # vs 
  sig_marg^2
  1/(1-rho^2)
  #
  
  
  
  d_tmp <- CJ( id = 1:5, t_obs = 1:3 )
  set.seed(1)
  d_tmp[, obs := rnorm(.N)]
  d_tmp[, y := 2 + 2 * obs + rnorm(.N)]
  d_tmp <- d_tmp[!(id %in% c(1) & t_obs == 3), ]
  
  
  
  d_tmp[, ix_tmp := 1:.N, keyby = id]
  d_ix <- d_tmp[, .(ix_s = .I[1], ix_e = .I[.N]), by=id]
  
  i <- 1
  mu <- rep(NA_real_, nrow(d_tmp))
  j <- 1
  k <- 2
  phi <- 0.5
  for(j in 1:length(unique(d_tmp$id))){
    N_id = d_ix[j, ix_e - ix_s + 1]
    cat("Number in grp ", j, " is ", N_id , "\n")
    
    cat("Set mu at index ", d_ix[j, ix_s] , "\n")
    mu[i] <-  d_tmp[d_ix[j, ix_s], 1 + 2 * age]
    
    epsi <- rep(NA, N_id)
    
    cat("Set first innovation for this group \n")
    epsi[1] = d_tmp[d_ix[j, ix_s], y] - mu[i]
    
    for(k in 2:N_id){
      cat("Calc innovation ",k," for group ", j, "\n")
      epsi[k] = d_tmp[d_ix[j, ix_s + k - 1], y] - mu[]
        
      i <- i + 1
      
      d_tmp[i, mu := d_tmp[i + j - 1, mu] + phi*epsi[k-1]]
      
    }
  }
  d_tmp[]
  
 
  
  library(brms)
  data("LakeHuron")
  LakeHuron <- as.data.frame(LakeHuron)
  brms::make_stancode(x ~ ar(p = 2), data = LakeHuron)
  brms::make_standata(x ~ ar(p = 2), data = LakeHuron)
  
  fit <- brm(x ~ ar(p = 2), data = LakeHuron)
  summary(fit)
  
  
   
  
  library(data.table)
  set.seed(1)
  d <- CJ(
    trt = 1:3,
    time = c(0, 0.25, 0.5, 0.75, 1)
  )
  d[, trt := factor(trt)]
  d[, time := factor(time)]
  d[, mu_tru := rnorm(.N)]
  d[, y := rnorm(.N, mu_tru, 1)]
  
  model.matrix(y ~ trt * time, data = d)
  
}

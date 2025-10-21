
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
  if(is.null(l_spec$tfu)){ tfu <- 1 } else { tfu <- l_spec$tfu }
  
  trt_opts <- 1:length(l_spec$p_trt_alloc)
  
  N_ic <- length(is:ie)
  # Progression from pre to post assuming soc occurs for everyone
  d <- CJ(
    # interim id
    ic = ic, 
    # unit id
    id = is:ie,
    # spirometry observation point
    t_obs = l_spec$t_sprty_obs
    )
  
  # age distribution - truncated log normal - restrict age to (age_lwr,age_upr)
  p_lt_age_lwr <- plnorm(
    l_spec$age_lwr, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig, lower.tail = T)
  p_gt_age_upr <- plnorm(
    l_spec$age_upr, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig, lower.tail = F)
  # sample u in [p0, 1)
  u <- p_lt_age_lwr + runif(N_ic, 0, (1 - (p_lt_age_lwr + p_gt_age_upr))) 
  baseline_age <- qlnorm(u, meanlog = l_spec$age_mu, sdlog = l_spec$age_sig)
  d[, age0 := baseline_age[id]]
  
  # trt allocation
  trt_alloc <- trt_opts[rep(trt_opts[as.logical(l_spec$p_trt_alloc)], len = N_ic)]
  d[, trt := trt_alloc[id]]
  
  # mean
  d[, mu :=
      l_spec$b_0 +
      l_spec$b_ln_age * log(age0) +
      l_spec$b_decl * t_obs +
      l_spec$b_trt[trt]
    ]     
  
  # Covariance matrix for AR(1) structure
  Sigma <- outer(
    1:length(l_spec$t_sprty_obs), 
    1:length(l_spec$t_sprty_obs), 
    function(i, j){
      l_spec$Sigma_s^2 * l_spec$Sigma_rho^abs(i - j)
  } )
  
  # outcome
  i <- 600
  for(i in 1:N_ic){
    
    # d[id == i, y := as.numeric(
    #   rmvt(n = 1, sigma = Sigma, df = l_spec$mvt_df, delta = d[id == i, mu]))]
    
    d[id == i, y := as.numeric(
      rmvnorm(n = 1, mean = d[id == i, mu], sigma = Sigma))]
  }
  
  # missing
  d[, y_mis := rbinom(.N, 1, l_spec$pr_ymis)]
  
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
  l_spec$b_decl <- -1
  # coefficient on log-age
  l_spec$b_ln_age <- -14
  # l_spec$b_ln_age <- -2
  # treatment effects
  l_spec$b_trt <- c(0, 0, 0) 
  
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
  l_spec$t_sprty_obs <- c(0, 0.25, 0.5, 0.75, 1)
  
  # refine 
  # number of spirometry events per year (mixture of poisson)
  # so that the number of observations is unconstrained to specific visits
  # which is probably more likely to be the reality
  # plot_spirom_data <- function(r = 2, p = 0.9){
  #   x <- 0:25
  #   y = dnbinom(x, size = r, prob = p)
  #   plot(x, y, type = "h", xlim=c(0, 25), xaxt='n' )
  #   axis(side = 1, at = x )
  # }
  # plot_spirom_data(r = 1.2, p = 0.3)
  # number of evts (negbin but will truncate to force N > 0) 
  # l_spec$num_sprty_obs_nb_r <- 1.2 
  # l_spec$num_sprty_obs_nb_p <- 0.3
  
  # perits - ignore.
  
  # pt level variation in baseline
  # l_spec$u_s0 <- 5
  # l_spec$u_s1 <- 1
  # # correlation between pre and post
  # l_spec$u_rho <- -0.1
  
  # MVN seems to restrictive for these data. Adopt MVT AR1
  
  # Degrees of freedom for multivariate t 
  # as mvt_df -> Inf mvt -> mvn
  # (3-10) heavier tails
  # l_spec$mvt_df <- 8
  
  # AR(1) correlation structure
  l_spec$Sigma_rho <- 0.85      # correlation between adjacent time points
  l_spec$Sigma_s <- 6       # standard deviation at each time point
  

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

prototype_1_data <- function(){
  
  l_spec <- get_prototype_cfg()
  str(l_spec)

  d <- get_sim03_trial_data(l_spec)

  # plot_1(d)
  
  
  # assume treatment is randomized at baseline to individuals 
  # sampled from the same population
  # run a standard pre/post (ancova) analysis
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  
  m1 <- cmdstanr::cmdstan_model("stan/sim03-v06.stan")  
  
  
  d_mod <- dcast(d, id + age0 + trt ~ t_obs, value.var = "y")
  ld <- list(
    N = d_mod[, .N],
    J = length(l_spec$t_sprty_obs),
    y = d_mod[, 4:8],
    age = d_mod$age0,
    t_obs = l_spec$t_sprty_obs,
    trt = d_mod$trt
  )
  foutname_1 <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"),
    "-sim-", 1)
  
  f_1 <- m1$sample(
    ld, iter_warmup = 1000, iter_sampling = 1000,
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
  f_1$summary(
    variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho"))
  
  # alternatively could use:
  d_mod2 <- copy(d)
  d_mod2[, trt := factor(trt)]
  f_1_f <- nlme::gls(y ~ log(age0) + trt + t_obs, d_mod2,
             correlation = nlme::corAR1(form = ~ 1 | id))
  summary(f_1_f)
  
  # f_1_optim <- m1$optimize(data = ld, jacobian = TRUE)
  # f_1_a <- m1$laplace(data = ld, mode = f_1_optim, draws = 2000)
  # f_1_a$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho", "nu"))

  # f_1_b <- m1$variational(data = ld)
  # f_1_b$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho", "nu"))
  # 
  # f_1_c <- m1$pathfinder(data = ld)
  # f_1_c$output()
  # f_1_c$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho"))
  
  
  
  d_post_1 <- data.table(f_1$draws(
    variables = c("b0", "b_ln_age", "b_t", "b_trt_raw"), 
    format = "matrix"))
  
  d_fig_1 <- melt(d_post_1, measure.vars = names(d_post_1))
  
  ggplot(d_fig_1, aes(x = value)) +
    geom_density() +
    ggh4x::facet_wrap2(~variable, ncol = 1, scales = "free_x")
}


plot_1 <- function(d){
  
  # d[, .(mu = mean(y)), keyby = .(trt, t_obs)]
  # min_y <- min(floor(d$y))
  # max_y <- max(ceiling(d$y))
  
  d_fig <- copy(d)
  d_fig[, trt := factor(trt)]
  
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
  p2
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
  
  
   
}

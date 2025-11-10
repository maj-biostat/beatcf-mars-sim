
if(exists("prefix_r")){
  source(paste0(prefix_r, "/libs.R"))
  source(paste0(prefix_r, "/data.R"))
} else {
  source("./R/libs.R")
  source("./R/data.R")
}


get_sim04_trial_data <- function(
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
  
  d_i[, mu0 := l_spec$b_0 + l_spec$b_la * log(age0) ]
  # what we would actually observe and therefore model
  d_i[, y0 := rnorm(.N, mu0, l_spec$sigma_mu0)]
  
  d <- CJ(
    # unit id
    id = 1:N_ic,
    # spirometry observation point - excludes baseline
    t_id = l_spec$t_sprty_obs
  )
  d <- d[d_i, on = "id"]
  
  # Build basis
  # Should generalise this but haven't got around to it. 
  # Trying to make the basis such that it is easy to enter the treatment
  # effects. At the moment it works for observation at every month 1:12 inclusive
  # with four segments. So, for example, we have a linear slope that covers month
  # 1 to 3 inclusive, then from 3 to 6, 6 to 9 and 9 to 12. 
  n_pw_seg <- nrow(l_spec$basis_ix_seg)
  k <- 1
  for (k in 1:n_pw_seg) {
    # within-interval linear ramp scaled to [0,1]
    start <- l_spec$basis_ix_seg[k, 1]
    end   <- l_spec$basis_ix_seg[k, 2]
    seg_length <- end - start
    
    d[, paste0("B", k) :=
        fifelse(t_id < start, 0,
                fifelse(t_id >= end, 1,
                        (t_id - start) / seg_length))]
  }
  
  
  #--- Compute linear predictor (mean FEV1) ---
  d[, T2 := fifelse(trt == 2, 1, 0)]
  d[, T3 := fifelse(trt == 3, 1, 0)]
  d[, mu :=
      mu0 +
      (l_spec$b_t[1] * B1 + l_spec$b_t[2] * B2 + l_spec$b_t[3] * B3 + l_spec$b_t[4] * B4) +
      (l_spec$b_T[1] * T2 + l_spec$b_T[2] * T3) +
      (l_spec$b_tT[1,1] * B1*T2 + l_spec$b_tT[1,2] * B2*T2 + l_spec$b_tT[1,3] * B3*T2 + l_spec$b_tT[1,4] * B4*T2) +
      (l_spec$b_tT[2,1] * B1*T3 + l_spec$b_tT[2,2] * B2*T3 + l_spec$b_tT[2,3] * B3*T3 + l_spec$b_tT[2,4] * B4*T3)
  ]
  
  
  
  
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
  
  d[, ia := NA_integer_]
  d[, t_anlys := NA_real_]
  d[, id := rep(is:ie, each = l_spec$n_sprty_obs)]
  d[, t_fu := t0 + (t_id/l_spec$n_sprty_obs)*365]
  
  setcolorder(
    d,
    c("ic", "id", "t_id", "t0", "t_fu", "age0", 
      "trt", "mu0", "y0"))
  
  d
}




get_sim04_prototype_cfg <- function(){
  
  f_cfg_test <- file.path("./etc/sim04/cfg-sim04-v01.yml")
  cfg_test <- config::get(file = f_cfg_test)
  
  l_spec <- list()
  
  l_spec$desc <- "Test"
  l_spec$nsim <- 1
  l_spec$nex <- 1
  
  # N by analysis
  l_spec$N <- sum(cfg_test$N_pt)
  
  l_spec$pt_per_day <-   cfg_test$pt_per_day
  l_spec$ramp_up_days <-  cfg_test$ramp_up_days
  
  # number of spirometry followup (excl baseline)
  l_spec$n_sprty_obs <- cfg_test$n_sprty_obs
  l_spec$t_sprty_obs <- 1:l_spec$n_sprty_obs
  # assume linearity within intervals (over integer follow up visits)
  # basis segments - each row indexs the start/stop followup visit
  l_spec$basis_ix_seg <- matrix(
    unlist(cfg_test$basis_ix_seg),
    ncol = 2, byrow = T)
  
  
  l_spec$age_mu  <- cfg_test$age_mu
  l_spec$age_sig <- cfg_test$age_sig
  l_spec$age_lwr <- cfg_test$age_lwr
  l_spec$age_upr <- cfg_test$age_upr
  
  # trt alloc
  l_spec$p_trt_alloc <- unlist(cfg_test$trt)/length(unlist(cfg_test$trt))
  
  # parameters
  l_spec$b_0 <- cfg_test$b_0
  l_spec$b_la <- cfg_test$b_la
  l_spec$b_t <- unlist(cfg_test$b_t)
  l_spec$b_t <- sapply(1:length(l_spec$b_t), function(i) eval(parse(text = l_spec$b_t[i])))
  l_spec$b_T <- unlist(cfg_test$b_T)
  l_spec$b_tT <- matrix(
    unlist(cfg_test$b_tT),
    ncol = length(l_spec$b_t), byrow = T
  )
  
  # baseline measure
  l_spec$sigma_mu0 <- cfg_test$sigma_mu0
  
  # resid
  l_spec$sigma <- cfg_test$sigma
  l_spec$rho <- cfg_test$rho
  
  
  
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
  
  
  l_spec
}



prototype_1_data <- function(){
  
  l_spec <- get_sim04_prototype_cfg()
  
  str(l_spec)

  d <- get_sim04_trial_data(l_spec)

  d_fig <- d[, .(mu = mean(y)), keyby = .(trt, t_id)]
  ggplot(d_fig, aes(x = t_id, y = mu, col = factor(trt))) +
    geom_line() +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 70:90) +
    theme(legend.position = "bottom")
#
  d_fig <- copy(d)
  d_fig[, age0 := round(age0)]
  d_fig[, age := age0 + (t_id/l_spec$n_sprty_obs) ]
  d_fig[id == 1]
  
  ggplot(d_fig, aes(x = age, y = y)) +
    geom_line(aes(group = id), lwd = 0.2) +
    theme(legend.position = "bottom") +
    facet_wrap(~trt, ncol = 1)
  #
  
  # same (uncorrelated random effects) via lme
  f_4 <- nlme::gls(
    y ~ y0  + log(age0) + 
      B1 + B2 + B3 + B4 +
      T2 + T3 +
      B1:T2 + B2:T2 + B3:T2 + B4:T2 +
      B1:T3 + B2:T3 + B3:T3 + B4:T3,
    correlation = nlme::corAR1(form = ~ 1 | id),
    data = d); summary(f_4)
  
  
  
  l_spec <- get_sim04_prototype_cfg()
  r <- pbapply::pblapply(
    X=1:100, cl = 4, FUN=function(ix) {
      
      d <- get_sim04_trial_data(l_spec)
      
      f_4 <- nlme::gls(
        y ~ y0  + log(age0) + 
          B1 + B2 + B3 + B4 +
          T2 + T3 +
          B1:T2 + B2:T2 + B3:T2 + B4:T2 +
          B1:T3 + B2:T3 + B3:T3 + B4:T3,
        correlation = nlme::corAR1(form = ~ 1 | id),
        data = d); 
      
      s <- summary(f_4)
      
      # extract residual sd and autocorrel
      aux <- exp(as.numeric(f_4$modelStruct$corStruct))
      rho <- (aux - 1) / (aux + 1)
      
      res <- c(
        coef(f_4),
        sigma = s$sigma,
        rho = rho
      )
      res
    }
  )
  
  res <- do.call(rbind, r)
  
  d_tbl <- data.table(
    par = colnames(res),
    mu = round(colMeans(res), 3)
  )
  d_tbl
  
  
  l_spec <- get_sim04_prototype_cfg()
  d <- get_sim04_trial_data(l_spec)
  f_4 <- nlme::gls(
    y ~ y0  + log(age0) + 
      # time trend
      B1 + B2 + B3 + B4 +
      # trt
      T2 + T3 +
      # trt x time
      B1:T2 + B2:T2 + B3:T2 + B4:T2 +
      B1:T3 + B2:T3 + B3:T3 + B4:T3,
    correlation = nlme::corAR1(form = ~ 1 | id),
    data = d); 
  
  mu <- coef(f_4)
  S <- vcov(f_4)
  
  d_post <- data.table(
    rmvnorm(1e5, mu, sigma = S)
  )
  
  d_post_trt <- d_post[, .(
    delta_2_1 = `B1:T2` + `B2:T2` + `B3:T2` + `B4:T2`,
    delta_3_1 = `B1:T3` + `B2:T3` + `B3:T3` + `B4:T3`)]
  
  d_fig <- melt(d_post_trt, measure.vars = names(d_post_trt))
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_density() +
    facet_wrap(~variable, scales = "free")
  #
  
  
  # d_post_4 <- data.table(
  #   rmvnorm(1e5, fixef(f_4), sigma = vcov(f_4))
  # )
  # names(d_post_4) <- paste0("b_", 0:4)
  # d_post_4_long <- melt(d_post_4, measure.vars = names(d_post_4))
  # 
  # 
  # ggplot(d_post_4_long, aes(x = value, group = type, col = type)) +
  #   geom_density() +
  #   facet_wrap(~variable, scales = "free")
  #
  
  
  
  # X <- model.matrix(~ t_obs + log(age0) + trt, data = d)
  
  
  # brms::make_stancode(
  #   y ~ t_obs + log(age0) + trt + (1 | id), data = d, family = "normal"
  # )
  
  
  # output_dir_mcmc <- paste0(getwd(), "/tmp")
  
  # m1 <- cmdstanr::cmdstan_model("stan/sim04-v01.stan")  
  
  
  # all those that have completed followup to 1 year
  # lsd <- list(
  #   N = d[, .N],
  #   N_S = d[, .GRP, by = id][, .N],
  #   id = d$id,
  #   t = d$t_obs,
  #   P = ncol(X) - 1,
  #   X = X[, -1],
  #   y = d[, y]
  # )
  # 
  # foutname_1 <- paste0(
  #   format(Sys.time(), format = "%Y%m%d%H%M%S"),
  #   "-sim-", 1)
  # 
  # f_1 <- m1$sample(
  #   lsd, iter_warmup = 1000, iter_sampling = 1000,
  #   parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = T,
  #   max_treedepth = 10,
  #   output_dir = output_dir_mcmc,
  #   output_basename = foutname_1
  # )
  # 
  # f_1$summary(variables = c("b_0", "b", "s_u0", "s_u1", "s_e"))
  
  
  # m2 <- cmdstanr::cmdstan_model("stan/sim04-v02.stan")  
  # lsd <- brms::make_standata(
  #   y ~ t_obs + log(age0) + trt + (1 | id), data = d, family = "normal"
  # )
  # f_2 <- m2$sample(
  #   lsd, iter_warmup = 1000, iter_sampling = 1000,
  #   parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = T,
  #   max_treedepth = 10,
  #   output_dir = output_dir_mcmc,
  #   output_basename = foutname_1
  # )
  # f_2$summary(variables = c("b_Intercept", "b", "sd_1", "sigma"))
  
  # m3 <- brms::make_stancode(
  #   y ~ t_obs + log(age0) + trt + (1| id) + (0 + t_obs | id), data = d, family = "normal"
  # )
  # m3 <- cmdstanr::cmdstan_model("stan/sim04-v03.stan")
  # lsd <- brms::make_standata(
  #   y ~ t_obs + log(age0) + trt + (1| id) + (0 + t_obs | id), data = d, family = "normal"
  # )
  # f_3 <- m3$sample(
  #   lsd, iter_warmup = 1000, iter_sampling = 1000,
  #   parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = T,
  #   max_treedepth = 10,
  #   output_dir = output_dir_mcmc,
  #   output_basename = foutname_1
  # )
  # f_3$summary(variables = c("b_Intercept", "b", "sd_1", "sd_2", "sigma"))
  
  # same (uncorrelated random effects) via lme
  # f_4 <- nlme::lme(y ~ t_obs + log(age0) + trt,
  #                  random = list(id = pdDiag(~ 1 + t_obs)),
  #                  data = d)
  # 
  # d_post_1 <- data.table(
  #   f_1$draws(variables = c("b_0", "b"), format = "matrix")
  # )
  # names(d_post_1) <- paste0("b_", 0:4)
  # d_post_4 <- data.table(
  #   rmvnorm(1e5, fixef(f_4), sigma = vcov(f_4))
  # )
  # names(d_post_4) <- paste0("b_", 0:4)
  # d_post_1_long <- melt(d_post_1, measure.vars = names(d_post_1))
  # d_post_4_long <- melt(d_post_4, measure.vars = names(d_post_4))
  # 
  # d_fig <- rbind(
  #   cbind(type = "stan", d_post_1_long),
  #   cbind(type = "lme", d_post_4_long)
  # )
  # ggplot(d_fig, aes(x = value, group = type, col = type)) +
  #   geom_density() +
  #   facet_wrap(~variable, scales = "free")
  #
  
  
  
  # f_1_optim <- m1$optimize(data = ld, jacobian = TRUE)
  # f_1_a <- m1$laplace(data = ld, mode = f_1_optim, draws = 2000)
  # f_1_a$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho", "nu"))

  # f_1_b <- m1$variational(data = ld)
  # f_1_b$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho", "nu"))
  # 
  # f_1_c <- m1$pathfinder(data = ld)
  # f_1_c$output()
  # f_1_c$summary(variables = c("b0", "b_ln_age", "b_t", "b_trt_raw", "sigma", "rho"))
  
  
  
  # d_post_1 <- data.table(f_1$draws(variables = c("b0", "b"),  format = "matrix"))
  # 
  # d_fig_1 <- melt(d_post_1, measure.vars = names(d_post_1))
  # 
  # ggplot(d_fig_1, aes(x = value)) +
  #   geom_density() +
  #   ggh4x::facet_wrap2(~variable, ncol = 2, scales = "free_x")
}




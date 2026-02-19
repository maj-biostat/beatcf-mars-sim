# Experiment with independent set of models with reduced linear predictor.

source("./R/init.R")
source("./R/data.R")
source("./R/data-sim12.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim12"
  args[2] = "./sim12/cfg-sim12-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# Log setup
f_log_sim <- file.path("./logs", "log-sim.txt")
log_appender(appender_file(f_log_sim))
log_info("*** START UP ***")
# log_threshold(TRACE)

f_cfgsc <- file.path("./etc", args[2])
l_spec <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(l_spec))

ix <- 1


output_dir_mcmc <- paste0(getwd(), "/tmp")

# # maps frequentist to bayesian model parameter names
# par_map <- paste0("b_", 0:13)
# names(par_map) <- c(
#   "(Intercept)", "y0", "log(age0)", "t_obs0.5", "t_obs0.75", 
#   "t_obs1", "trt2", "trt3", "t_obs0.5:trt2", "t_obs0.75:trt2", 
#   "t_obs1:trt2", "t_obs0.5:trt3", "t_obs0.75:trt3", "t_obs1:trt3"
# )



# Main trial loop.
run_trial <- function(
    ix,
    l_spec,
    return_posterior = F
){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # Get enrolment times for arbitrary large set of patients
  # Simpler to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  
  # events per day
  lambda = l_spec$pt_per_day
  # ramp up over x months 
  rho = function(t) pmin(t/l_spec$ramp_up_days, 1)
  
  # lambda = 0.57
  # # ramp up over x months 
  # rho = function(t) pmin(t/120, 1)
  # 
  # rr <- unlist(pblapply(1:100, cl = 4, FUN=function(ii){
  #   ttt <- get_enrol_time(sum(l_spec$N_pt), lambda, rho)
  #   max(ttt)
  # }))
  # mean(rr) / 365
  
  # day of enrolment
  loc_t0 <- get_enrol_time(sum(l_spec$N_pt), lambda, rho)
  
  # d_fig <- data.table(id = 1:length(loc_t0), loc_t0, loc_t)
  # d_fig <- melt(d_fig, id.vars = "id")
  # 
  # ggplot(d_fig, aes(x = value, y = id, group = variable)) +
  #   geom_step()

  
  
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ic <- 1 # interim number
  N_analys <- length(l_spec$N_pt)
  
  
  # # posterior summaries
  # g_par <- c("b0", "y0", "log(age0)", "B1", "B2", "B3", "B4", "T2", 
  #            "T3", "B1:T2", "B2:T2", "B3:T2", "B4:T2", "B1:T3", "B2:T3", "B3:T3", 
  #            "B4:T3", "delta_2_1", "delta_3_1")
  # g_par_delta <- c("delta_2_1", "delta_3_1")
  #   
  # d_post_smry_1 <- CJ(
  #   ic = 1:N_analys,
  #   par = factor(c(g_par), levels = c(g_par))
  # )
  # d_post_smry_1[, mu := NA_real_]
  # d_post_smry_1[, med := NA_real_]
  # d_post_smry_1[, se := NA_real_]
  # d_post_smry_1[, q_025 := NA_real_]
  # d_post_smry_1[, q_975 := NA_real_]
  
  # decisions 
  # g_rule_type <- c("ni",  "fut")
  # 
  # d_pr_dec <- CJ(
  #   ic = 1:N_analys,
  #   rule = factor(g_rule_type),
  #   par = factor(g_par_delta),
  #   p = NA_real_,
  #   dec = NA_integer_
  # )
  
  # store all simulated trial pt data
  d_all <- data.table()
  # 
  # 
  # if(return_posterior){
  #   d_post_all <- data.table()
  # }
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " cohort ", l_spec$ic)
    
    # next chunk of data on pts.
    if(l_spec$ic == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N_pt[l_spec$ic] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N_pt[l_spec$ic] - 1
    }
    
    # id and time
    l_spec$t0 <- loc_t0[l_spec$is:l_spec$ie]
    
    # We are assuming that the analysis takes place on pt having reached endpoint
    
    d <- get_sim12_cohort(l_spec)
    
    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    if(l_spec$ie == sum(l_spec$N_pt)){
      log_info("Trial ", ix, " final analysis, using all pt")
      t_now <- d_all[, max(t_fu)]
      
      l_mod <- get_sim12_stan_data(d_all)
      
    } else {
      # t0 is the entry time and is repeated 
      t_now <- d_all[, max(t0)]
      
      # include all those who have reached the first follow up at the time of 
      # the interim analysis - we don't need to wait for them to complete follow
      # up in order to use them in the analysis - they still give us data on the 
      # progression of fev1 up to their latest visit. mfit also uses this 
      # perspective basically including pt if they have been on trt for at 
      # least 4 weeks.
      incl_ids <- d_all[t0 + 365 < t_now, id]
      
      l_mod <- get_sim12_stan_data(d_all[id %in% incl_ids])
      
      
    }
    
    # see help on msm2Surv and (jackson paper on multistate to see data structure
    # requirements)
    # d_tmp <- melt(d_all[, .(id, age, start, stop, trt, state)], id.vars = c("id", "age", "trt", "state"))
    # d_tmp <- d_tmp[variable == "start"]
    # # state needs to be numeric/factor
    # d_tmp[, state := factor(state)]
    # permitted transitions (diag ignored)
    # Q <- rbind(c(1, 1), 
    #            c(1, 1))
    # 
    # msm::msm2Surv(
    #   data = d_tmp, 
    #   subject = "id", 
    #   time = "value", 
    #   state= "state",
    #   covs = c("trt", "age"),
    #   Q=Q)
    
    f1 <- joint(
      formSurv =
        list(
          inla.surv(time = gap_time, event = evt_he) ~ ppfev_std + (1 | id),
          inla.surv(time = gap_time, event = evt_eh) ~ ppfev_std + trt + (1 | id)
        ),
      basRisk= "weibullsurv", id = "id", dataSurv = d_mod, control= list(config=TRUE)
    )
    
    d_smry <- rbind(
      cbind(
        model = "S1",
        parameter = rownames(summary(f1)$SurvEff[[1]]),
        data.table(summary(f1)$SurvEff[[1]], hr = T)
      ),
      cbind(
        model = "S1",
        parameter = rownames(summary(f1)$ReffListS[[1]]),
        data.table(summary(f1)$ReffListS[[1]])
      ),
      cbind(
        model = "S2",
        parameter = rownames(summary(f1)$SurvEff[[2]]),
        data.table(summary(f1)$SurvEff[[2]], hr = T)
      ),
      cbind(
        model = "S2",
        parameter = rownames(summary(f1)$ReffListS[[2]]),
        data.table(summary(f1)$ReffListS[[2]])
      ),
      fill = T
    )
    
    
    l_draws <- inla.posterior.sample(1e3, f1)
    l_draws[[1]]$hyperpar
    l_draws[[2]]$hyperpar
    l_draws[[3]]$hyperpar
    l_draws[[4]]$hyperpar
    
    # state transition probabilities
    
    # difference in probability of recovery by day 14 
    
    # predict for entire cohort probability of recovery at day 14 assuming assigned to each trt group
    # compute mean and difference between mean probability of recovery.
    
    # for B draws from the posterior, compute prob of recov at day 14 for each pt then take expectation
    # over the B draws this gives us a distribution for each group
    
    
    # draw from multivariate normal on all parameters
    d_post <- data.table(
      rmvnorm(1e4, mu, sigma = S)
    )
    # isolate treatment effects of interest for convenience
    d_post[, `:=`(
      delta_2_1 = `B1:T2` + `B2:T2` + `B3:T2` + `B4:T2`,
      delta_3_1 = `B1:T3` + `B2:T3` + `B3:T3` + `B4:T3`
      )]
    
    log_info("Trial ", ix, " fitted models ", l_spec$ic)
    
    if(return_posterior){
      d_post_all <- rbind(
        d_post_all,
        cbind(ic = l_spec$ic, d_post)
      )
    }
    
    
    d_post_long <- melt(d_post, measure.vars = names(d_post), variable.name = "par")
    
    
    # merge posterior summaries for current interim
    d_post_smry_1[
      d_post_long[, .(ic = l_spec$ic,
                      mu = mean(value),
                      med = median(value),
                      se = sd(value),
                      q_025 = quantile(value, prob = 0.025),
                      q_975 = quantile(value, prob = 0.975)
                      ), keyby = par], 
      on = .(ic, par), `:=`(
        mu = i.mu,
        med = i.med,
        se = i.se,
        q_025 = i.q_025,
        q_975 = i.q_975
        )]

    log_info("Trial ", ix, " extracted posterior ", l_spec$ic)
    
    
    # compute and merge the current probability and decision trigger status
    
    # The trial setup is such that if the intervention were beneficial, it would 
    # reduce the temporal decline of fev1. 
    
    # If an intervention arm is non-inferior then the fev at 12 months (landmark)
    # will decline by no more than the non-inferiority margin relative to the ctl 
    
    # If the probability that the fev on the intervention arm is worse (lower) 
    # than the non-inferiority margin is very high, then we will conclude 
    # inferiority
    
    d_pr_dec[
      
      rbind(
        d_post_long[par %in% g_par_delta, .(
          ic = l_spec$ic,
          rule = factor("ni", levels = g_rule_type),
          p = mean(value > l_spec$delta$ni),
          dec = as.integer(mean(value > l_spec$delta$ni) > l_spec$thresh$ni)
        ), keyby = par],
        d_post_long[par %in% g_par_delta, .(
          ic = l_spec$ic,
          rule = factor("fut", levels = g_rule_type),
          p = mean(value < l_spec$delta$fut),
          dec = as.integer(mean(value < l_spec$delta$fut) > l_spec$thresh$fut)
        ), keyby = par]
      ),
      
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec  
      )
    ]
    
    
    # For trial stopping we only consider the comparisons to the
    # soc. In order to stop the study both need to have been resolved (either 
    # NI or futility have been concluded).
    d_stop <- d_pr_dec[
      ic <= l_spec$ic & par %in% c(g_par_delta),
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]
    
    
    # Update allocation probabilities if decisions have been reached.
    # Subsequent enrolments get redirected to the remaining arms.
    if(all(d_stop$resolved)){
      log_info("Stop trial as all questions addressed ", ix)
      stop_enrol <- T
    } else if(any(d_stop$resolved)){

      # 2 vs 1 at timepoint 4
      if(d_stop[par == "delta_2_1", resolved]) {
        l_spec$p_trt_alloc[2] <- 0
        l_spec$p_trt_alloc <- l_spec$p_trt_alloc / sum(l_spec$p_trt_alloc)
      }

      # 3 vs 1 at timepoint 4
      if(d_stop[par == "delta_3_1", resolved]) {
        l_spec$p_trt_alloc[3] <- 0
        l_spec$p_trt_alloc <- l_spec$p_trt_alloc / sum(l_spec$p_trt_alloc)
      }
    }
    log_info("Trial ", ix, " allocation after cohort ", l_spec$ic, " alloc: ", paste0(l_spec$p_trt_alloc, collapse = ", "))
    
    
    # next interim
    l_spec$ic <- l_spec$ic + 1
    
    if(l_spec$ic > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- l_spec$ic - 1
  
  
  l_ret <- list(
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    d_pr_dec = d_pr_dec,
    stop_at = stop_at
  )
  
  if(return_posterior){
    l_ret$d_post_all <- copy(d_post_all)
  }
  # 
  
  return(l_ret)
}






run_sim12 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 3
  }
  
  # trt alloc - balanced over number of trts
  l_spec$b_trt <- unlist(l_spec$b_trt)
  l_spec$n_trt <- length(l_spec$b_trt)
  
  
  if(l_spec$nex > 0){
    log_info("Creating ", l_spec$nex, " example trials with full posterior")
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$nsim, size = l_spec$nex, replace = F))
  }
  return_posterior <- F
  str(l_spec)
  e = NULL
  ix <- 1
  
  ## LOOP -------
  
  
  log_info("Start simulation")
  
  r <- pbapply::pblapply(
    X=1:l_spec$nsim, cl = l_spec$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);
      
      if(ix %in% l_spec$ex_trial_ix){
        return_posterior = T  
      } else {
        return_posterior = F
      }
      
      ll <- tryCatch({
        run_trial(
          ix,
          l_spec,
          return_posterior = return_posterior
        )
      },
      error=function(e) {
        log_info("ERROR in MCLAPPLY LOOP (see terminal output):")
        message(" ERROR in MCLAPPLY LOOP " , e);
        log_info("Traceback (see terminal output):")
        message(traceback())
        stop(paste0("Stopping with error ", e))
      })
      
      # ll$d_all[, sum(N)]
      ll
    })
  
  
  log_info("Length of result set ", length(r))
  log_info("Sleep for 5 before processing")
  Sys.sleep(5)
  
  for(i in 1:length(r)){
    log_info("Element at index ",i, " is class ", class(r[[i]]))
    if(any(class(r[[i]]) %like% "try-error")){
      log_info("Element at index ",i, " has content ", r[[i]])  
    }
    log_info("Element at index ",i, " has names ", 
             paste0(names(r[[i]]), collapse = ", "))
  }
  
  
  d_pr_dec <- data.table()
  for(i in 1:length(r)){
    
    log_info("Appending d_pr_dec for result ", i)
    
    if(is.recursive(r[[i]])){
      d_pr_dec <- rbind(
        d_pr_dec,
        cbind(
          sim = i, r[[i]]$d_pr_dec
        ) 
      )  
    } else {
      log_info("Value for r at this index is not recursive ", i)
      message("r[[i]] ", r[[i]])
      message(traceback())
      stop(paste0("Stopping due to non-recursive element "))
    }
    
  }
  
  d_post_smry_1 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_1
  } ), idcol = "sim")
  
  
  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_all
  } ), idcol = "sim")
  
  
  d_post_all <- data.table(do.call(rbind, lapply(1:length(r), function(i){
    # if the sim contains full posterior (for example trial) then return
    if(!is.null(r[[i]]$d_post_all)){
      cbind(sim = i, r[[i]]$d_post_all)
    }
    
  } )))
  
  l <- list(
    cfg = l_spec,
    d_pr_dec = d_pr_dec, 
    d_post_smry_1 = d_post_smry_1,
    d_all = d_all,
    d_post_all = d_post_all
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim12/sim12-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim12 <- function(){
  log_info("run_none_sim12: Nothing doing here bud.")
}

main_sim12 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim12()












# f_2 <- gamm(y ~ y0_std + log(age0) + t_id * trt,
#             random = list(id = ~ 1),
#             correlation = corAR1(form = ~ as.numeric(t_id) | id),
#             data = d_mod, method = "REML")

# summary(f_2$gam)
# summary(f_2$lme)   

# d_new <- data.table(t_id = factor(12), 
#                     trt = factor(1), 
#                     y0_std = mean(d_mod$y0_std), 
#                     age0 = mean(d_mod$age0))
# pred_1 <- predict(f_2$gam, d_new, se.fit = TRUE)
# d_new$trt <- factor(2)
# pred_2 <- predict(f_2$gam, d_new, se.fit = TRUE)
# d_new$trt <- factor(3)
# pred_3 <- predict(f_2$gam, d_new, se.fit = TRUE)
# 
# delta_2_1_est <- pred_2$fit - pred_1$fit
# delta_2_1_se  <- sqrt(pred_1$se.fit^2 + pred_2$se.fit^2)
# 
# delta_3_1_est <- pred_3$fit - pred_1$fit
# delta_3_1_se  <- sqrt(pred_1$se.fit^2 + pred_3$se.fit^2)
# 
# d_post_long[, .(mu = mean(value)), keyby = par]
# c(delta_2_1_est, delta_3_1_est)




# d_fig_1 <- melt(d_post, measure.vars = names(d_post))
# 
# d_fig_1[, .(mean(value), sd(value)), keyby = variable]
# ggplot(d_fig_1, aes(x = value)) +
#   geom_density() +
#   facet_wrap(~variable, scales = "free_x")
# 
# d_fig_2 <- data.table(f_1$draws(
#   variables = c(
#     
#     "b_trt[2]", "b_trt[3]", "b_pre", "se"
#     
#   ),   
#   format = "matrix"))
# d_fig_2 <- melt(d_fig_2, measure.vars = names(d_fig_2))
# d_fig_2[, .(mean(value), sd(value)), keyby = variable]
# ggplot(d_fig_2, aes(x = value)) +
#   geom_density() +
#   facet_wrap(~variable, scales = "free_x")
#

# d_fig_3 <- rbind(
#   d_fig_1[variable %in% c("delta_2_1", "delta_3_1")]  ,
#   d_fig_2[variable %in% c("b_trt[2]", "b_trt[3]")]  
# )
# d_fig_3[variable %in% c("delta_2_1", "b_trt[2]"), par := "delta_2_1"]
# d_fig_3[variable %in% c("delta_3_1", "b_trt[3]"), par := "delta_3_1"]
# 
# ggplot(d_fig_3, aes(x = value, group = variable, col = variable)) +
#   geom_density() +
#   facet_wrap(variable~par, scales = "free_x")

# d_mu <- data.table(f_1$draws(
#   variables = c(
# 
#     c("mu")
# 
#   ),
#   format = "matrix"))
# d_mu <- melt(d_mu, measure.vars = names(d_mu))

# d_fig_1 <- copy(d_mod)
# d_fig_1[, y_pred := d_mu[, mean(value), keyby = variable][, V1]]
# ggplot(d_fig_1, aes(x = y_post, y = y_pred)) +
#   geom_point()



# snk <- capture.output(
#   f_1_mode <- m1$optimize(data = lsd$ld, jacobian = TRUE)
# )
# snk <- capture.output(
#   f_1 <- m1$laplace(data = lsd$ld)
# )


# f_1 <- m1$sample(
#   lsd$ld, iter_warmup = 1000, iter_sampling = 1000,
#   parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = T,
#   max_treedepth = 11,
#   output_dir = output_dir_mcmc,
#   output_basename = foutname
# )




# d_fig_1 <- melt(d_post, measure.vars = names(d_post))
# d_fig_1[, .(mean(value), sd(value)), keyby = variable]
# ggplot(d_fig_1, aes(x = value)) +
#   geom_density() +
#   facet_wrap(~variable, scales = "free_x")

# d_fig_2 <- data.table(f_1$draws(
#   variables = c("b", "se"),
#   format = "matrix"))
# d_fig_2 <- melt(d_fig_2, measure.vars = names(d_fig_2))
# d_fig_2[, .(mean(value), sd(value)), keyby = variable]
# ggplot(d_fig_2, aes(x = value)) +
#   geom_density() +
#   facet_wrap(~variable, scales = "free")

# b_0 <- as.numeric(f_1$draws(variables = c("b_0"), format = "matrix"))
# b <- f_1$draws(variables = c("b"), format = "matrix")
# X <- model.matrix(~ y_pre + factor(trt), data = d_mod)[, -1]

# d_fig_3 <- data.table(b_0 + t(X %*% t(b)))
# d_fig_3 <- data.table(
#   mu = apply(d_fig_3, 2, mean),
#   q_lwr = apply(d_fig_3, 2, FUN = quantile, prob = 0.025),
#   q_upr = apply(d_fig_3, 2, FUN = quantile, prob = 0.975)
# )

# d_fig_4 <- data.table(
#   mu_f = predict(ff)
# )

# d_fig_5 <- copy(d_mod[y_mis == 0])
# d_fig_5 <- cbind(d_fig_5, d_fig_3[, .(mu, q_lwr, q_upr)])
# d_fig_5 <- cbind(d_fig_5, d_fig_4)

# ggplot(d_fig_5, aes(x = id, y = y_post)) +
#   geom_point(size = 0.5) +
#   geom_point(aes(x = id, y = mu), col = 2, size = 1.5) +
#   geom_linerange(
#     aes(x = id, ymin = q_lwr, ymax = q_upr), , col = 2
#   ) +
#   geom_point(aes(x = id, y = mu_f), col = 3, size = 0.5) 

# ggplot(d_fig_5, aes(x = y_pre, y = y_post, col = factor(trt), group = trt)) +
#   geom_point(size = 0.5) +
#   geom_point(aes(y = mu), size = 1.5) 
# geom_linerange(
#   aes(x = id, ymin = q_lwr, ymax = q_upr), , col = 2
# ) +
# geom_point(aes(x = id, y = mu_f), col = 3, size = 0.5) 
#

# d_fig_2[, .(mean(value), sd(value)), keyby = variable]
# ggplot(d_fig_2, aes(x = value)) +
#   geom_density() +
#   facet_wrap(~variable, scales = "free_x")



source("R/libs.R")
source("R/init.R")
source("R/util.R")
source("R/data-sim15.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim15"
  args[2] = "./sim15/cfg-sim15-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}


f_log <- file.path("./logs", "log.txt")
log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")


# Log setup
f_log_sim <- file.path("./logs", "log-sim.txt")
log_appender(appender_file(f_log_sim))
log_info("*** START UP ***")
# log_threshold(TRACE)

f_cfgsc <- file.path("./etc", args[2])
l_spec <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(l_spec))


ix <- 1
m1 <- cmdstanr::cmdstan_model("stan/sim15-v02.stan")

output_dir_mcmc <- paste0(getwd(), "/tmp")



#' Main trial loop, sequentially creating patient sample, assigning trt and 
#' running analyses and decision processes.
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
  
  # day of enrolment
  loc_t0 <- sim_ipp_thinning(sum(l_spec$N_pt), lambda, rho)
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ic <- 1 # interim number
  N_analys <- length(l_spec$N_pt)
  
  # posterior summary for parameters from models
  d_post_smry_1 <- CJ(
    ic = 1:N_analys,
    par = l_spec$par_names
  )
  
  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, lo := NA_real_]
  d_post_smry_1[, hi := NA_real_]
  
  # trt effects
  d_trt_effects <- CJ(
    ic = 1:N_analys,
    par = c(
      "soc", "def", "dis", "delta_def", "delta_dis"
    ),
    mu = NA_real_,
    lwr = NA_real_,
    upr = NA_real_
  )
  
  # decisions 
  d_pr_dec <- CJ(
    ic = 1:N_analys,
    rule = c("ni", "fut"),
    trt = l_spec$trt_lab[2:3],
    p = NA_real_,
    dec = NA_integer_
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
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
    
    d <- get_sim15_cohort(l_spec)
    
    # averge duration of state 
    # d[, .N, by = .(id, state)][, mean(N), by = state]

    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # segment by state by adding run ids
    d_all[state == "H", bin := l_spec$v_lu_he_bin[day_in_state + 1L] ]
    d_all[state == "E", bin := l_spec$v_lu_eh_bin[day_in_state + 1L] ]
    d_all[, rlgrp := rleid(id, state, trt, bin)]
    
    # are we at the final analysis or interim?
    if(l_spec$ie == sum(l_spec$N_pt)){
      
      log_info("Trial ", ix, " final analysis, using all pt")
      t_now <- d_all[, max(t0 + l_spec$followup)]
      l_mod <- get_sim15_stan_data(d_all, l_spec)
      
    } else {
      # t0 is the entry time (note that this is repeated for each id if they 
      # have more than one event over the 12 months)
      t_now <- d_all[, max(t0)]
      
      # no longer restriction on pt having had one exacerbation, we include
      # all data
      incl_ids <- d_all[t0 + day_of_fu < t_now, unique(id)]
      l_mod <- get_sim15_stan_data(dd = copy(d_all[id %in% incl_ids]), l_spec)
      
    }
    
    f_1 <- m1$sample(
      l_mod, 
      iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
      parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain, 
      refresh = 100, show_exceptions = T,
      max_treedepth = 11
    )
    # f_1$summary(variables = l_spec$par_names_pre)
    
    d_post <- data.table(
      f_1$draws(
        format = "matrix",
        variables = l_spec$par_names_pre
      )
    )
    
    names(d_post) <- l_spec$par_names
    
    # posterior summary (eh model)
    d_post_smry_1[
      data.table(
        ic = l_spec$ic, 
        par = names(d_post),
        mu = colMeans(d_post),
        lo = apply(d_post, 2, function(z){
          quantile(z, prob = 0.025)
        }),
        hi = apply(d_post, 2, function(z){
          quantile(z, prob = 0.975)
        })
        ),
      on = .(ic, par), `:=`(
        mu = i.mu, lo = i.lo, hi = i.hi
      )
    ]
    
    # B_max = 300
    # N_pt = l_mod$N_id
    # ppfev_0 = unique(d_all$ppfev_0)
    d_res <- calc_trt_effect(
      d_post, B_max = 200, N_pt = l_mod$N_id, 
      # unique ppfev0 from sample
      ppfev_0 = unique(d_all$ppfev_0), l_spec)
    
    # d_fig <- melt(d_res, measure.vars = names(d_res))
    # 
    # ggplot(
    #   d_fig[variable %in% l_spec$trt_lab], 
    #        aes(x = value, group = variable)) +
    #   geom_density() + 
    #   facet_wrap(~variable, ncol = 1)
    
    # update treatment effects 
    d_trt_effects[
      data.table(
        ic = l_spec$ic, 
        par = names(d_res),
        mu = colMeans(d_res),
        lwr = apply(d_res, 2, function(z){quantile(z, prob = 0.025)}),
        upr = apply(d_res, 2, function(z){quantile(z, prob = 0.975)})
      ),
      on = .(ic, par), `:=`(
        mu = i.mu, lwr = i.lwr, upr = i.upr
      )
    ]
    
    # evaluate decision rules, namely does the rmst indicate a longer duration of 
    # recovery in the intervention group relative to the soc group that is 
    # above the level that we are willing to tolerate
    
    # NI
    d_res_def <- data.table(
      ic = l_spec$ic,  rule = "ni", trt = "def", 
      # hoping that any increase in recovery time is below the NI margin
      p = mean(d_res$delta_def < l_spec$dec_eh_delta_ni) 
    )
    # for NI decision, probability must exceed our evidential threshold
    d_res_def[, dec := as.integer(p > l_spec$dec_thresh_ni)]
    
    # same but for early discontinue
    d_res_dis <- data.table(
      ic = l_spec$ic,  rule = "ni", trt = "dis", 
      p = mean(d_res$delta_dis < l_spec$dec_eh_delta_ni)
    )
    d_res_dis[, dec := as.integer(p > l_spec$dec_thresh_ni)]
    
    d_pr_dec[
      rbind(
        d_res_def,
        d_res_dis
      ),
      on = .(ic, rule, trt), `:=`(
        p = i.p, dec = i.dec
      )
    ]
    
    
    # futility
    d_res_def <- data.table(
      ic = l_spec$ic,  rule = "fut", 
      trt = "def", 
      # unfortunately the revised schedule might be substantially greater than NI
      # in the case of exacerbation
      p = mean(d_res$delta_def > l_spec$dec_eh_delta_ni)
    )
    # for futility decision, probability must exceed our evidential thresholds
    d_res_def[, dec := as.integer(p > l_spec$dec_thresh_fut)]
    
    d_res_dis <- data.table(
      ic = l_spec$ic,  rule = "fut", 
      trt = "dis", 
      p = mean(d_res$delta_dis > l_spec$dec_eh_delta_ni) 
    )
    d_res_dis[, dec := as.integer(p > l_spec$dec_thresh_fut)]
    
    d_pr_dec[
      rbind(
        d_res_def,
        d_res_dis
      ),
      on = .(ic, rule, trt), `:=`(
        p = i.p, dec = i.dec
      )
    ]
    
    # evaluate status - have we stopped?
    d_stop <- d_pr_dec[
      ic <= l_spec$ic & trt %in% l_spec$trt_lab[2:3],
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(trt)]
    
    if(all(d_stop$resolved)){
      log_info("Stop trial as all questions addressed ", ix)
      stop_enrol <- T
    } else if(any(d_stop$resolved)){
      
      if(d_stop[trt == "def", resolved]) {
        l_spec$trt_active["def"] <- FALSE
      }
      
      if(d_stop[trt == "dis", resolved]) {
        l_spec$trt_active["dis"] <- FALSE
      }
    }
    
    log_info("Trial ", ix, " allocation after cohort ", 
             l_spec$ic, " alloc: ", paste0(l_spec$trt_active, collapse = ", "))
    
    
    # next interim
    l_spec$ic <- l_spec$ic + 1
    d_all[, `:=`(
      bin = NULL, rlgrp = NULL
    )]
    
    if(l_spec$ic > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- l_spec$ic - 1
  
  
  l_ret <- list(
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    # d_post_smry_2 = d_post_smry_2,
    # d_post_smry_3 = d_post_smry_3,
    d_trt_effects = d_trt_effects, 
    d_pr_dec = d_pr_dec,
    stop_at = stop_at,
    l_spec = l_spec
  )
  
  # 
  
  return(l_ret)
}





#' Entry point for running trial simulation in parallel.
run_sim15 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 3
  }
  
  
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
  
  l_spec$v_lu_he_bin <- l_spec$d_lu_he_bin$ix_bin
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
        run_trial(ix, l_spec, return_posterior = return_posterior )
      },
      error=function(e) {
        log_info("ERROR in MCLAPPLY LOOP (see terminal output):")
        message(" ERROR in MCLAPPLY LOOP " , e); message(traceback())
        stop(paste0("Stopping with error ", e))
      })
      
      # d_tbl <- ll$d_all[, .N, by = .(id, state)]
      # all_ids <- unique(ll$d_all$id)
      # all_states <- c("H", "E")
      # 
      # d_grid <- CJ(id = all_ids, state = all_states)
      # d_tbl_full <- d_tbl[d_grid, on = .(id, state)]
      # d_tbl_full[is.na(N), N := 0]
      # # if a pt has one exacerbation but is censored (ie we never see the 
      # # recovery time because it exceeds fu) then the below summary will show
      # # a min of 1 for both E and H
      # table(d_tbl_full$state, d_tbl_full$N)
      # # naive average duration of E
      # ll$d_all[state == "E", .(.N, dur_mu = mean(dur)), keyby = trt]
      # # units entering into each analysis (enrolment defer)
      # ll$d_all[t0 %in% ll$l_spec$t0_last, .SD[1], keyby = id]
      # # posterior
      # ll$d_post_smry_1[ic == ll$stop_at, .(par, mu = round(mu, 3), lo = round(lo, 3), hi = round(hi, 3))]
      # # rmst
      # ll$d_post_smry_2[ic == ll$stop_at, .(trt, rmst_mu = round(rmst_mu, 2), lo = round(rmst_lo, 2), hi = round(rmst_hi, 2))]
      # # trt effect, decision
      # ll$d_pr_dec[ic <= ll$stop_at]
      #
      
      
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
  
  # d_post_smry_2 <- rbindlist(lapply(1:length(r), function(i){ 
  #   r[[i]]$d_post_smry_2
  # } ), idcol = "sim")
  # 
  # d_post_smry_3 <- rbindlist(lapply(1:length(r), function(i){ 
  #   r[[i]]$d_post_smry_3
  # } ), idcol = "sim")
  
  
  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_all
  } ), idcol = "sim")
  
  d_analys_sets <- rbindlist(lapply(1:length(r), function(i){ 
    data.table(t0_last = r[[i]]$l_spec$t0_last)
  } ), idcol = "sim")
  
  
  l <- list(
    l_spec = l_spec,
    d_pr_dec = d_pr_dec, 
    d_post_smry_1 = d_post_smry_1,
    # d_post_smry_2 = d_post_smry_2,
    # d_post_smry_3 = d_post_smry_3,
    d_all = d_all,
    d_analys_sets = d_analys_sets
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim15/sim15-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  message("fname is ", fname)
  
  log_info("Saving results file ", fname)
  
  message("saving... ", fname)
  
  qs::qsave(l, file = fname)
  
  
  message("saved")
}

run_none_sim15 <- function(){
  log_info("run_none_sim15: Nothing doing here bud.")
}

main_sim15 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim15()












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



source("R/libs.R")
source("R/init.R")
source("R/util.R")
source("R/data-sim16.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim16"
  args[2] = "./sim16/cfg-sim16-v01.yml"
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
m1 <- cmdstanr::cmdstan_model("stan/sim16-v01.stan")

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
  
  # t_600 <- numeric(100)
  # for(i in 1:100){
  #   t_600[i] = sim_ipp_thinning(sum(l_spec$N_pt), lambda, rho)[600]
  # }
  # hist(t_600)
  
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
  
  # decisions 
  d_pr_dec <- CJ(
    ic = 1:N_analys,
    rule = c("ni", "fut"),
    trt = l_spec$trt_lab[2:3],
    p = NA_real_,
    dec = NA_integer_
  )
  
  # enrolment 
  d_enrol <- data.table()
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  d_pri_par <- data.table()
  d_post_par <- data.table()
  
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
    
    d <- get_sim16_cohort(l_spec)
    
    # averge duration of state 
    # d[, .N, by = .(id, state)][, mean(N), by = state]

    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # are we at the final analysis or interim?
    if(l_spec$ie == sum(l_spec$N_pt)){
      
      log_info("Trial ", ix, " final analysis, using all pt")
      t_now <- d_all[, max(t0 + l_spec$followup + l_spec$decision_fu)]
      # only those with events contribute
      l_mod <- get_sim16_stan_data(dd = copy(d_all[evt == 1]), l_spec)
      
      d_enrol <- rbind(
        d_enrol,
        d_all[evt == 1, .SD[.N], keyby = id][
          , .(ic = l_spec$ic, t_now, id)]
      )
      
    } else {
      # t0 is the entry time (note that this is repeated for each id if they 
      # have more than one event over the 12 months)
      t_now <- d_all[, max(t0)]
      
      # we include participants that have had an exacerbation and for whom
      # we have completed the subsequent 90 day follow up
      incl_ids <- d_all[evt == 1 & t0 + day_of_exac + l_spec$decision_fu < t_now, unique(id)]
      l_mod <- get_sim16_stan_data(dd = copy(d_all[id %in% incl_ids]), l_spec)
    
      d_enrol <- rbind(
        d_enrol,
        d_all[evt == 1 & t0 + day_of_exac + l_spec$decision_fu < t_now, .SD[.N], keyby = id][
          , .(ic = l_spec$ic, t_now, id)]
      )
      
    }
    
    f_1 <- m1$sample(
      l_mod, 
      iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
      parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain, 
      refresh = 0, show_exceptions = T,
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
    
    # kableExtra::kable(d_post_smry_1[ic == l_spec$ic], format = "simple", digits = 2)
    
    
    if(return_posterior){
      
      # compute priors
      l_mod$prior_only <- 1
      f_1_pri <- m1$sample(
        l_mod, 
        iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
        parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain, 
        refresh = 0, show_exceptions = T,
        max_treedepth = 11
      )
      d_pri <- data.table(
        f_1_pri$draws(
          format = "matrix",
          variables = l_spec$par_names_pre
        )
      )
      names(d_pri) <- l_spec$par_names
      
      d_pri_par <- rbind(
        d_pri_par,
        cbind(ic = l_spec$ic, d_pri)
      )
      
      # posterior
      d_post_par <- rbind(
        d_post_par,
        cbind(ic = l_spec$ic, d_post)
      )
      
    }
    
    
    # evaluate decision rules, making decisions on the basis of expected 
    # period in exacerbation state 
    
    # NI
    d_res_def <- data.table(
      ic = l_spec$ic,  rule = "ni", trt = "def", 
      # hoping that any increase in recovery time is below the NI margin
      p = mean(d_post$b_def < l_spec$dec_eh_delta_ni) 
    )
    # for NI decision, probability must exceed our evidential threshold
    d_res_def[, dec := as.integer(p > l_spec$dec_thresh_ni)]
    
    # same but for early discontinue
    d_res_dis <- data.table(
      ic = l_spec$ic,  rule = "ni", trt = "dis", 
      p = mean(d_post$b_dis < l_spec$dec_eh_delta_ni)
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
      p = mean(d_post$b_def > l_spec$dec_eh_delta_ni)
    )
    # for futility decision, probability must exceed our evidential thresholds
    d_res_def[, dec := as.integer(p > l_spec$dec_thresh_fut)]
    
    d_res_dis <- data.table(
      ic = l_spec$ic,  rule = "fut", 
      trt = "dis", 
      p = mean(d_post$b_dis > l_spec$dec_eh_delta_ni) 
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
        
        log_info("Trial ", ix, 
                 " def stopped no further enrolment into this arm, analy id ", 
                 l_spec$ic)
      }
      
      if(d_stop[trt == "dis", resolved]) {
        l_spec$trt_active["dis"] <- FALSE
        
        log_info("Trial ", ix, 
                 " dis stopped no further enrolment into this arm, analy id ", 
                 l_spec$ic)
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
  
  
  # lobstr::obj_size(d_w)
  # lobstr::obj_size(d_all)
  
  l_ret <- list(
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    d_pr_dec = d_pr_dec,
    d_pri_par = d_pri_par,
    d_post_par = d_post_par,
    d_enrol = d_enrol,
    stop_at = stop_at,
    l_spec = l_spec
  )
  
  # 
  
  return(l_ret)
}





#' Entry point for running trial simulation in parallel.
run_sim16 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 3
  }
  
  l_spec <- cfg_update(l_spec)
  
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
      
      ll
    })
  
  
  Sys.sleep(10)
  
  # l_ret <- list(
  #   d_all = d_all,
  #   d_post_smry_1 = d_post_smry_1,
  #   # d_post_smry_2 = d_post_smry_2,
  #   # d_post_smry_3 = d_post_smry_3,
  #   d_trt_effects = d_trt_effects, 
  #   d_pr_dec = d_pr_dec,
  #   stop_at = stop_at,
  #   l_spec = l_spec
  # )
  # 
  
  # data from each simulated trial
  d_w <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_w
  } ), idcol = "sim")
  
  d_post_smry_1 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_1
  } ), idcol = "sim")
  
  d_trt_effects <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_trt_effects
  } ), idcol = "sim")
  
  d_pr_dec <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_pr_dec
  } ), idcol = "sim")
  
  d_pri_par <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_pri_par
  } ), idcol = "sim")
  
  d_pri_res <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_pri_res
  } ), idcol = "sim")
  
  d_post_par <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_par
  } ), idcol = "sim")
  
  d_post_res <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_res
  } ), idcol = "sim")
  
  d_enrol <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_enrol
  } ), idcol = "sim")
  
  l <- list(
    l_spec = l_spec,
    d_w = d_w,
    d_post_smry_1 = d_post_smry_1,
    d_trt_effects = d_trt_effects,
    d_pr_dec = d_pr_dec,
    d_pri_par = d_pri_par,
    d_pri_res = d_pri_res,
    d_post_par = d_post_par,
    d_post_res = d_post_res,
    d_enrol = d_enrol
  )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim16/sim16-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  message("fname is ", fname)
  qs::qsave(l, file = fname)
  
  
  message("saved")
}

run_none_sim16 <- function(){
  log_info("run_none_sim16: Nothing doing here bud.")
}

main_sim16 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim16()







source("R/libs.R")
source("R/init.R")
source("R/util.R")
source("R/data-sim18.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim18"
  args[2] = "./sim18/cfg-sim18-v01.yml"
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
m_1 <- cmdstanr::cmdstan_model("stan/sim18-v01.stan")

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
  l_spec$t_0 <- sim_ipp_thinning(sum(l_spec$N_pt), lambda, rho)
  
  # t_600 <- numeric(100)
  # for(i in 1:100){
  #   t_600[i] = sim_ipp_thinning(sum(l_spec$N_pt), lambda, rho)[600]
  # }
  # hist(t_600)
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ic <- 1 # interim number
  N_analys <- length(l_spec$N_pt)
  
  # summary for model parameters
  d_post_smry_1 <- CJ(
    ic = 1:N_analys,
    par = l_spec$full_pars
  )

  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, lo := NA_real_]
  d_post_smry_1[, hi := NA_real_]
  
  # summary sop
  d_post_smry_2 <- CJ(
    ic = 1:N_analys,
    state = l_spec$state_lab,
    trt = l_spec$trt_lab,
    day = l_spec$visit_days[-1]
  )
  d_post_smry_2[, mu := NA_real_]
  d_post_smry_2[, lo := NA_real_]
  d_post_smry_2[, hi := NA_real_]
  
  # summary time in state, difference in time in state
  d_post_smry_3 <- CJ(
    ic = 1:N_analys,
    state = l_spec$state_lab,
    par = c(l_spec$trt_lab, l_spec$delta_lab)
  )
  d_post_smry_3[, mu := NA_real_]
  d_post_smry_3[, lo := NA_real_]
  d_post_smry_3[, hi := NA_real_]
  
  # decisions 
  d_post_smry_4 <- CJ(
    ic = 1:N_analys,
    rule = c("ni", "fut"),
    par = l_spec$delta_lab,
    p = NA_real_,
    dec = NA_integer_
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  # full posterior 
  d_post_par <- data.table()
  d_post_sop <- data.table()
  
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
    
    
    # We are assuming that the analysis takes place on pt having reached endpoint
    
    d <- sim18_cohort(l_spec)$d_obs

    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # are we at the final analysis or interim?
    if(l_spec$ie == sum(l_spec$N_pt)){
      
      log_info("Trial ", ix, " final analysis, using all pt")
      # only those with events contribute
      l_mod <- sim18_stan_data(dd = copy(d_all), l_spec)
      
    } else {
      # t_0 is the entry time (note that this is repeated for each id if they 
      # have more than one event over the 12 months)
      t_now <- d_all[, max(t_0)]
      
      # we include participants that have had an exacerbation and for whom
      # we have completed the subsequent 90 day follow up
      incl_ids <- d_all[t_0 + l_spec$followup_dec < t_now, unique(id)]
      l_mod <- sim18_stan_data(dd = copy(d_all[id %in% incl_ids]), l_spec)
    
    }
    
    f_1_optim <- m_1$optimize(data = l_mod, jacobian = TRUE)
    f_1 <- m_1$laplace(data = l_mod, mode = f_1_optim, draws = 2000, refresh = 0)
    
    # f_1 <- m_1$sample(
    #   l_mod,
    #   iter_warmup = l_spec$mcmc_warmup, iter_sampling = l_spec$mcmc_iter,
    #   parallel_chains = l_spec$mcmc_chain, chains = l_spec$mcmc_chain,
    #   refresh = 0, show_exceptions = T,
    #   max_treedepth = 11
    # )
    
    m_post <- f_1$draws(variables = l_spec$smry_pars, format = "matrix")
    
    l_spec_mod <- copy(l_spec)
    d_sop <- rbindlist(lapply(1:nrow(m_post), function(ii){

      # calculate sop based on parameter estimates at their post mean
      l_spec_mod$alpha <- as.numeric(m_post[ii , c("a[1]", "a[2]")])
      l_spec_mod$b_trt <- as.numeric(m_post[ii , c("b_trt[1]", "b_trt[2]", "b_trt[3]")])  
      l_spec_mod$b_prev <- as.numeric(m_post[ii , c("b_prev[1]", "b_prev[2]", "b_prev[3]")])  
      l_spec_mod$b_time_1 <- as.numeric(m_post[ii , c("b_time_1")])  
      l_spec_mod$b_time_2 <- as.numeric(m_post[ii , c("b_time_2")]) 
      l_spec_mod$b_gap <- as.numeric(m_post[ii , c("b_gap[1]", "b_gap[2]")])
      l_spec_mod$b_trt_time <- as.numeric(m_post[ii , c("b_trt_time[1]", "b_trt_time[2]", "b_trt_time[3]")])  

      names(l_spec_mod$b_trt) <- l_spec_mod$trt_lab
      names(l_spec_mod$b_trt_time) <- l_spec_mod$trt_lab

      d_tmp <- sim18_sop(days = l_spec$visit_days[-1], l_spec_mod)
      d_tmp
      
    }), idcol = "id_draw" )
    
    
    if(return_posterior){
      
      # posterior
      d_post_par <- rbind(
        d_post_par,
        cbind(ic = l_spec$ic, m_post)
      )
      
      # posterior
      d_post_sop <- rbind(
        d_post_sop,
        cbind(ic = l_spec$ic, d_sop)
      )
      
    }
    
    
    d_sop[, gap := day - data.table::shift(day, type = "lag"), keyby = .(id_draw, trt)]

    d_sop[, `:=`(
      prop_none = none * gap,
      prop_mild = mild * gap,
      prop_severe = severe * gap)]
    
    
    # days in each state
    d_dur <- melt(
      d_sop[day > 0, .(
      none = sum(none * gap),
      mild = sum(mild * gap),
      severe = sum(severe * gap)
      ), keyby = .(trt, id_draw)], id.vars = c("trt", "id_draw"), 
      variable.name = "state", value.name = "days")
    
    d_dur <- dcast(d_dur, state + id_draw ~ trt, value.var = "days")
    d_dur[, delta_def := def - soc]
    d_dur[, delta_dis := dis - soc]
    
    # sd(d_dur$delta_def)
    # sd(d_dur$delta_dis)
    
    # posterior summary (did the model recover the data generating pars)
    d_post <- data.table(m_post)
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
    
    
    # posterior summary - what are the estimated sops
    d_post_smry_2[
      melt(d_sop[day > 0, .(id_draw, trt, day, none, mild, severe)], 
           id.vars = c("id_draw", "trt", "day"), variable.name = "state")[
             , .(ic = l_spec$ic,
                 mu = mean(value), 
                 lo = quantile(value, prob = 0.025),
                 hi = quantile(value, prob = 0.975)),
             keyby = .(state, trt, day)
           ],
      on = .(ic, state, trt, day), `:=`(
        mu = i.mu, lo = i.lo, hi = i.hi
      )
    ]
    
    # posterior summary - days in state and treatment effects
    d_post_smry_3[
      melt(d_dur, id.vars = c("state", "id_draw"), variable.name = "par")[
        , .(
          ic = l_spec$ic,
          mu = mean(value),
          lo = quantile(value, prob  = 0.025),
          hi = quantile(value, prob  = 0.975)
        ), keyby = .(state, par)
      ],
      on = .(ic, par, state), `:=`(
        mu = i.mu, lo = i.lo, hi = i.hi
      )
    ]
    
    
    # posterior summary - evaluating decision rules on the differences
    d_post_smry_4[
      melt(d_dur[
        state == "none", .(id_draw, delta_def, delta_dis)], 
        id.vars = c("id_draw"), variable.name = "par")[
        , .(
          ic = l_spec$ic,
          rule = "ni",
          # probability that difference in number of symptom free days is above 
          # some nominal reduction
          p = mean(value > -l_spec$dec_delta_ni) ,
          dec = NA_real_
        ), keyby = .(par)
      ],
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec
      )
    ]
    d_post_smry_4[ic == l_spec$ic & rule == "ni", dec := as.integer(p > l_spec$dec_thresh_ni)]
    
    d_post_smry_4[
      melt(d_dur[state == "none", .(id_draw, delta_def, delta_dis)], id.vars = c("id_draw"), variable.name = "par")[
        , .(
          ic = l_spec$ic,
          rule = "fut",
          # probability that difference in number of symptom free days is above 
          # some nominal reduction
          p = mean(value < -l_spec$dec_delta_ni) ,
          dec = NA_real_
        ), keyby = .(par)
      ],
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec
      )
    ]
    d_post_smry_4[ic == l_spec$ic & rule == "fut", dec := as.integer(p > l_spec$dec_thresh_fut)]
    
    
    
    # evaluate status - have we stopped?
    d_stop <- d_post_smry_4[
      ic <= l_spec$ic,
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]
    
    if(all(d_stop$resolved)){
      log_info("Stop trial as all questions addressed ", ix)
      stop_enrol <- T
    } else if(any(d_stop$resolved)){
      
      if(d_stop[par == "delta_def", resolved]) {
        l_spec$trt_active["def"] <- FALSE
        
        log_info("Trial ", ix, 
                 " def stopped no further enrolment into this arm, analy id ", 
                 l_spec$ic)
      }
      
      if(d_stop[par == "delta_dis", resolved]) {
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
    d_post_smry_2 = d_post_smry_2,
    d_post_smry_3 = d_post_smry_3,
    d_post_smry_4 = d_post_smry_4,
    d_post_par = d_post_par,
    d_post_sop = d_post_sop,
    stop_at = stop_at,
    l_spec = l_spec
  )
  
  
  
  # 
  
  return(l_ret)
}





#' Entry point for running trial simulation in parallel.
run_sim18 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 3
  }
  
  l_spec <- update_sim18_cfg(l_spec)
  
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
  
  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_all
  } ), idcol = "sim")
  
  d_post_smry_1 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_1
  } ), idcol = "sim")
  
  d_post_smry_2 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_2
  } ), idcol = "sim")
  
  d_post_smry_3 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_3
  } ), idcol = "sim")
  
  d_post_smry_4 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_4
  } ), idcol = "sim")
  
  d_post_par <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_par
  } ), idcol = "sim")
  
  d_post_sop <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_sop
  } ), idcol = "sim")
  
  l <- list(
    l_spec = l_spec,
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    d_post_smry_3 = d_post_smry_3,
    d_post_smry_4 = d_post_smry_4,
    d_post_par = d_post_par,
    d_post_sop = d_post_sop
  )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim18/sim18-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  message("fname is ", fname)
  qs::qsave(l, file = fname)
  
  
  message("saved")
}

run_none_sim18 <- function(){
  log_info("run_none_sim18: Nothing doing here bud.")
}

main_sim18 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim18()





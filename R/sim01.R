# Experiment with independent set of models with reduced linear predictor.

source("./R/init.R")
source("./R/data-sim01.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim01"
  args[2] = "./sim01/cfg-sim01-v01.yml"
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
g_cfgsc <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(g_cfgsc))

ix <- 1
m1 <- cmdstanr::cmdstan_model("stan/sim01-v01.stan")

output_dir_mcmc <- paste0(getwd(), "/tmp")



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
  lambda = 1.0
  # ramp up over 12 months 
  rho = function(t) pmin(t/360, 1)
  
  loc_t0 <- get_sim01_enrol_time(sum(l_spec$N), lambda, rho)
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ia <- 1 # interim number
  N_analys <- length(l_spec$N)
  
  # posterior summaries
  g_par_p = c("p_1", "p_2", "p_3")
  g_par_rd = c("rd_2_1", "rd_3_1", "rd_3_2")
  d_post_smry_1 <- CJ(
    ia = 1:N_analys,
    par = factor(c(g_par_p, g_par_rd))
  )
  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, med := NA_real_]
  d_post_smry_1[, se := NA_real_]
  d_post_smry_1[, q_025 := NA_real_]
  d_post_smry_1[, q_975 := NA_real_]
  
  
  # decisions 
  g_rule_type <- c("ni",  "inf")
  
  d_pr_dec <- CJ(
    ia = 1:N_analys,
    rule = factor(g_rule_type),
    par = factor(g_par_rd),
    p = NA_real_,
    dec = NA_integer_
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  if(return_posterior){
    d_post_all <- data.table()
  }
  
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " analysis ", l_spec$ia)
    
    # next chunk of data on pts.
    if(l_spec$ia == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ia] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ia] - 1
    }
    
    # id and time
    l_spec$t0 = loc_t0[l_spec$is:l_spec$ie]
    
    # We are assuming that the analysis takes place on pt having reached endpoint
    
    d <- get_sim01_trial_data(l_spec)
    
    log_info("Trial ", ix, " new data generated ", l_spec$ia)
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_sim01_stan_data(d_all)
    
    lsd$ld$pri_a <- l_spec$prior$a
    lsd$ld$pri_b_trt <- l_spec$prior$btrt
    lsd$ld$pri_b_lung <- l_spec$prior$blung
    lsd$ld$pri_b_pseudo <- l_spec$prior$bpseudo
    
    # str(lsd$ld)
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", ix, "-intrm-", l_spec$ia)
    
    
    
    # fit model - does it matter that I continue to fit the model after the
    # decision is made...?
    snk <- capture.output(
      f_1 <- m1$sample(
        lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
        parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
        max_treedepth = 11,
        output_dir = output_dir_mcmc,
        output_basename = foutname
      )
    )
    
    log_info("Trial ", ix, " fitted models ", l_spec$ia)
    
    # extract posterior - marginal probability of outcome by trt group
    # higher values of p indicate higher risk of treatment failure
    d_post <- data.table(f_1$draws(
      variables = c(
        
        "p_1", "p_2", "p_3", "rd_2_1", "rd_3_1", "rd_3_2"
        
      ),   
      format = "matrix"))
    
    if(return_posterior){
      d_post_all <- rbind(
        d_post_all,
        cbind(ia = l_spec$ia, d_post)
      )
    }
    
    d_post_long <- melt(d_post, measure.vars = names(d_post), variable.name = "par")
    
    # d_fig <- copy(d_post_long)
    # d_fig[par %like% "p", par_type := "prob"]
    # d_fig[par %like% "rd", par_type := "rd"]
    # 
    # p1 <- ggplot(d_fig[par %like% "p", ],
    #              aes(x = value, group = par, col = par)) +
    #   geom_density()
    # p2 <- ggplot(d_fig[par %in% c("rd_2_1", "rd_3_1"), ],
    #              aes(x = value, group = par, col = par)) +
    #   geom_density()
    # p1/p2

    
    # merge posterior summaries for current interim
    d_post_smry_1[
      d_post_long[, .(ia = l_spec$ia,
                      mu = mean(value),
                      med = median(value),
                      se = sd(value),
                      q_025 = quantile(value, prob = 0.025),
                      q_975 = quantile(value, prob = 0.975)
                      ), keyby = par], 
      on = .(ia, par), `:=`(
        mu = i.mu,
        med = i.med,
        se = i.se,
        q_025 = i.q_025,
        q_975 = i.q_975
        )]

    log_info("Trial ", ix, " extracted posterior ", l_spec$ia)
    
    
    # compute and merge the current probability and decision trigger status
    
    # The trial setup is such that if the intervention were beneficial, it would 
    # reduce the event rate and we are therefore considering positive 
    # non-inferiority margins.
    # If an intervention arm is non-inferiority then the event rate on that arm 
    # will increase by no more than the non-inferiority margin.
    # If the probability that the event rate on the intervention arm is lower 
    # than the non-inferiority margin is very high, then we will concluded 
    # non-inferiority.
    # Contrastingly, for inferiority, if the probability that the risk 
    # difference is greater than the non-inferiority margin is high, then we 
    # will conclude inferiority. 
    
    d_pr_dec[
      
      rbind(
        d_post_long[par %in% g_par_rd, .(
          ia = l_spec$ia,
          rule = factor("ni", levels = g_rule_type),
          p = mean(value < l_spec$delta$ni),
          dec = as.integer(mean(value < l_spec$delta$ni) > l_spec$thresh$ni)
        ), keyby = par],
        d_post_long[par %in% g_par_rd, .(
          ia = l_spec$ia,
          rule = factor("inf", levels = g_rule_type),
          p = mean(value > l_spec$delta$inf),
          dec = as.integer(mean(value > l_spec$delta$inf) > l_spec$thresh$inf)
        ), keyby = par]
      ),
      
      on = .(ia, rule, par), `:=`(
        p = i.p, dec = i.dec  
      )
    ]
    
    log_info("Trial ", ix, " compared to thresholds ", l_spec$ia)
    
    # For trial stopping we only consider the comparisons to the
    # soc. In order to stop the study both need to have been resolved, either
    # NI or inferior.
    d_stop <- d_pr_dec[
      ia <= l_spec$ia & par %in% c("rd_2_1", "rd_3_1"), 
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]
    
    
    # Update allocation probabilities if decisions have been reached.
    # Subsequent enrolments get redirected to the remaining arms.
    if(any(d_stop$resolved)){
    
      if(d_stop[par == "rd_2_1", resolved]) {
        l_spec$p_trt_alloc[2] <- 0  
        l_spec$p_trt_alloc <- l_spec$p_trt_alloc / sum(l_spec$p_trt_alloc)
      } 
      
      if(d_stop[par == "rd_3_1", resolved]) {
        l_spec$p_trt_alloc[3] <- 0  
        l_spec$p_trt_alloc <- l_spec$p_trt_alloc / sum(l_spec$p_trt_alloc)
      } 
    }
    
    
    if(all(d_stop$resolved)){
      log_info("Stop trial as all questions addressed ", ix)
      stop_enrol <- T    
    }
    
    log_info("Trial ", ix, " updated allocation control ", l_spec$ia)
    
    # next interim
    l_spec$ia <- l_spec$ia + 1
    
    if(l_spec$ia > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- l_spec$ia - 1
  
  
  l_ret <- list(
    # data collected in the trial
    d_all = d_all[, .(y = sum(y), .N), keyby = .(ia, lung, pseudo, trt)],
    
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






run_sim01 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 3
  }
  
  l_spec <- list()
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  l_spec$p_lung_alloc <- g_cfgsc$lung_alloc
  l_spec$p_pseudo_alloc_lt70 <- g_cfgsc$pseudo_alloc_lt70
  l_spec$p_pseudo_alloc_ge70 <- g_cfgsc$pseudo_alloc_ge70
  
  # trt alloc
  l_spec$p_trt_alloc <- unlist(g_cfgsc$trt)/length(unlist(g_cfgsc$trt))
  
  
  # model parameters
  # model params (are on the risk scale)
  l_spec$bmu <- g_cfgsc$bmu
  l_spec$blung <- unlist(g_cfgsc$blung)
  l_spec$bpseudo <- unlist(g_cfgsc$bpseudo)
  l_spec$btrt <- unlist(g_cfgsc$btrt)
  
  # priors on log-odds/log-or
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$a <- unlist(g_cfgsc$pri_a)
  l_spec$prior$blung <- unlist(g_cfgsc$pri_blung)
  l_spec$prior$bpseudo <- unlist(g_cfgsc$pri_bpseudo)
  l_spec$prior$btrt <- unlist(g_cfgsc$pri_btrt)
  
  # decision hurdle
  l_spec$delta <- list()
  l_spec$delta$ni <- g_cfgsc$dec_delta_ni
  l_spec$delta$inf <- g_cfgsc$dec_delta_inf
  
  # evidentiary requirement
  l_spec$thresh <- list()
  l_spec$thresh$ni <- unlist(g_cfgsc$dec_thresh_ni)
  l_spec$thresh$inf <- unlist(g_cfgsc$dec_thresh_inf)
  
  l_spec$nex <- g_cfgsc$nex
  if(l_spec$nex > 0){
    log_info("Creating ", l_spec$nex, " example trials with full posterior")
    l_spec$ex_trial_ix <- sort(sample(1:g_cfgsc$nsim, size = l_spec$nex, replace = F))
  }
  return_posterior <- F
  str(l_spec)
  
  
  ## LOOP -------
  
  
  e = NULL
  log_info("Start simulation")
  r <- pbapply::pblapply(
    X=1:g_cfgsc$nsim, cl = g_cfgsc$mc_cores, FUN=function(ix) {
      # X=1:5, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
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
    cfg = g_cfgsc,
    d_pr_dec = d_pr_dec, 
    d_post_smry_1 = d_post_smry_1,
    d_all = d_all,
    d_post_all = d_post_all
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim01/sim01-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim01 <- function(){
  log_info("run_none_sim01: Nothing doing here bud.")
}

main_sim01 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim01()



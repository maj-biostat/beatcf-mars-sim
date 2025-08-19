# Experiment with independent set of models with reduced linear predictor.

source("./R/init.R")
# same data generating process
source("./R/data.R")
source("./R/data-sim01.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim00"
  args[2] = "./sim00/cfg-sim00-v06.yml"
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

output_dir_mcmc <- paste0(getwd(), "/tmp")

ix <- 1

# Main trial loop.
run_trial00 <- function(
    ix,
    l_spec,
    return_posterior = F
){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # loop controls
  N_total <- sum(l_spec$N)
  
  if(return_posterior){
    d_post_all <- data.table()
  }
  
  l_spec$is <- 1
  l_spec$ie <- N_total
  
  l_spec$N <- N_total
  
  d <- get_trial_data(l_spec)
  # d[]
  
  pri_a <- l_spec$prior$p[1]
  pri_b <- l_spec$prior$p[2]
  d_post <- data.table(
    p_1 = rbeta(2000, pri_a + d[trt == 1, sum(y)], pri_b + d[trt == 1, .N] - d[trt == 1, sum(y)]),
    p_2 = rbeta(2000, pri_a + d[trt == 2, sum(y)], pri_b + d[trt == 2, .N] - d[trt == 2, sum(y)]),
    p_3 = rbeta(2000, pri_a + d[trt == 3, sum(y)], pri_b + d[trt == 3, .N] - d[trt == 3, sum(y)])
  )
  d_post[, rd_2_1 := p_2 - p_1]
  d_post[, rd_3_1 := p_3 - p_1]
  d_post[, rd_3_2 := p_3 - p_2]
  d_post <- melt(d_post, measure.vars = names(d_post), variable.name = "par")
  
  # ggplot(d_post, aes(x = value)) +
  #   geom_density() +
  #   ggh4x::facet_wrap2(~variable) +
  #   theme_bw()
  #
  
  
  d_post_smry_1 <- d_post[, .(
    mu = mean(value),
    med = median(value),
    se = sd(value), 
    q_025 = quantile(value, prob = 0.025), 
    q_975 = quantile(value, prob = 0.975)
  ), keyby = .(par)]
  
  
  
  d_pr_dec <- rbind(
    d_post[par %like% "rd_", .(
      rule = "inf",
      p = mean(value > l_spec$delta$inf)
      ), keyby = .(par)],
    d_post[par %like% "rd_", .(
      rule = "ni",
      p = mean(value < l_spec$delta$ni)
    ), keyby = .(par)]
  )
  d_pr_dec[rule == "inf", dec := as.integer(p > l_spec$thresh$inf)]
  d_pr_dec[rule == "ni", dec := as.integer(p > l_spec$thresh$ni)]
  
  l_ret <- list(
    # data collected in the trial
    d_all = d[, .(y = sum(y), .N), keyby = .(ia, lung, pseudo, trt)],
    
    d_post_smry_1 = d_post_smry_1,
    
    d_pr_dec = d_pr_dec
  )
  
  if(return_posterior){
    l_ret$d_post_all <- copy(d_post)
  }
  # 
  
  l_ret
  
  
}

run_sim00 <- function(){
  
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
  l_spec$prior$p <- unlist(g_cfgsc$pri_p)
  
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
        run_trial00(
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
  
  d_pr_dec <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_pr_dec
  } ), idcol = "sim")
  
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
  
  fname <- paste0("data/sim00/sim00-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
  
  
  
}

run_none_sim00 <- function(){
  log_info("run_none_sim00: Nothing doing here bud.")
}

main_sim00 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim00()


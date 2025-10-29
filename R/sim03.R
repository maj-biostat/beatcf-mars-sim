# Experiment with independent set of models with reduced linear predictor.

source("./R/init.R")
source("./R/data.R")
source("./R/data-sim03.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim03"
  args[2] = "./sim03/cfg-sim03-v01.yml"
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
m1 <- cmdstanr::cmdstan_model("stan/sim03-v06b.stan")

output_dir_mcmc <- paste0(getwd(), "/tmp")

# maps frequentist to bayesian model parameter names
par_map <- paste0("b_", 0:13)
names(par_map) <- c(
  "(Intercept)", "y0", "log(age0)", "t_obs0.5", "t_obs0.75", 
  "t_obs1", "trt2", "trt3", "t_obs0.5:trt2", "t_obs0.75:trt2", 
  "t_obs1:trt2", "t_obs0.5:trt3", "t_obs0.75:trt3", "t_obs1:trt3"
)



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
  #   ttt <- get_enrol_time(sum(l_spec$N), lambda, rho)
  #   max(ttt)
  # }))
  # mean(rr) / 365
  
  # day of enrolment
  loc_t0 <- get_enrol_time(sum(l_spec$N), lambda, rho)
  
  # d_fig <- data.table(id = 1:length(loc_t0), loc_t0, loc_t)
  # d_fig <- melt(d_fig, id.vars = "id")
  # 
  # ggplot(d_fig, aes(x = value, y = id, group = variable)) +
  #   geom_step()
  #
  
  
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ic <- 1 # interim number
  N_analys <- length(l_spec$N)
  
  
  # posterior summaries
  g_par_beta = paste0("b_", 0:13)
  g_par_delta = paste0(rep(paste0("delta_", 2:3, "_1_"), each = 4), 1:4)
  
  g_par_delta_pri <- c("delta_2_1_4", "delta_3_1_4")
    
  d_post_smry_1 <- CJ(
    ic = 1:N_analys,
    par = factor(c(g_par_beta, g_par_delta))
  )
  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, med := NA_real_]
  d_post_smry_1[, se := NA_real_]
  d_post_smry_1[, q_025 := NA_real_]
  d_post_smry_1[, q_975 := NA_real_]
  
  # decisions 
  g_rule_type <- c("ni",  "inf")
  
  d_pr_dec <- CJ(
    ic = 1:N_analys,
    rule = factor(g_rule_type),
    par = factor(g_par_delta_pri),
    p = NA_real_,
    dec = NA_integer_
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  if(return_posterior){
    d_post_all <- data.table()
    d_freq <- data.table()
  }
  
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " cohort ", l_spec$ic)
    
    # next chunk of data on pts.
    if(l_spec$ic == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
    }
    
    # id and time
    l_spec$t0 <- loc_t0[l_spec$is:l_spec$ie]
    
    # We are assuming that the analysis takes place on pt having reached endpoint
    
    d <- get_sim03_trial_data(l_spec)
    
    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    if(l_spec$ie == sum(l_spec$N)){
      log_info("Trial ", ix, " final analysis, using all pt")
      t_now <- d_all[, max(t_fu)]
      d_mod <- copy(d_all)
      
      d_all[is.na(ia) & id %in% d_mod$id, ia := l_spec$ic]
      d_all[is.na(t_anlys) & id %in% d_mod$id, t_anlys := t_now]
      
    } else {
      # t0 is the entry time and is repeated 
      t_now <- d_all[, max(t0)]
      
      # include all those who have reached the first follow up at the time of 
      # the interim analysis
      incl_ids <- d_all[t_id == 1 & t_fu <= t_now, id]
      d_mod <- d_all[id %in% incl_ids]
      
      # set which analysis the cohort enters
      d_all[is.na(ia) & id %in% d_mod$id, ia := l_spec$ic]
      # set the analysis time for each unit
      d_all[is.na(t_anlys) & id %in% d_mod$id, t_anlys := t_now]
    }
    
    d_mod[, y0_std := scale(y0)]
    d_mod[, trt := factor(trt)]
    d_mod[, t_obs := factor(t_obs)]
    
    X <- model.matrix(~ y0_std + log(age0) + t_obs*trt, data = d_mod)
    
    # all those that have completed followup to 1 year
    lsd <- list(
      N = d_mod[, .N],
      P = ncol(X)-1,
      X = X[, -1],
      y = d_mod[, y],
      N_1 = d_mod[t_id == 1, .N],
      N_2 = d_mod[t_id == 2, .N],
      N_3 = d_mod[t_id == 3, .N],
      N_4 = d_mod[t_id == 4, .N],
      ix_1 = d_mod[t_id == 1, which = T],
      ix_2 = d_mod[t_id == 2, which = T],
      ix_3 = d_mod[t_id == 3, which = T],
      ix_4 = d_mod[t_id == 4, which = T]
    )
    
    # lsd$ld$pri_b_0 <- l_spec$prior$pri_b_0
    # lsd$ld$pri_b_mu <- l_spec$prior$pri_b_mu
    # lsd$ld$pri_b_sd <- l_spec$prior$pri_b_sd
    # lsd$ld$pri_se <- l_spec$prior$pri_b_se
    
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", ix, "-intrm-", l_spec$ic)
    
    snk <- capture.output(
      
      f_1 <- m1$sample(
        lsd, iter_warmup = 1000, iter_sampling = 1000,
        parallel_chains = 1, chains = 1, refresh = 100, show_exceptions = F,
        max_treedepth = 10,
        output_dir = output_dir_mcmc,
        output_basename = foutname,
        init = list(
          list(rho_un = 0, sigma = 1, b_0 = 100, b = rep(0, lsd$P))
        )
      )
      
    )
    
    log_info("Trial ", ix, " fitted models ", l_spec$ic)
    
    # extract posterior - marginal probability of outcome by trt group
    # higher values of p indicate higher risk of treatment failure
    m_post <- f_1$draws(variables = c("b"), format = "matrix")
    
    # treatment effects
    # | compare   | mnth  | contrast   |
    # |  2 vs 1   | 3     | b_6        |
    # |  2 vs 1   | 6     | b_6 + b_8  |
    # |  2 vs 1   | 9     | b_6 + b_9  |
    # |  2 vs 1   | 12    | b_6 + b_10 |
    # |  3 vs 1   | 3     | b_7        |
    # |  3 vs 1   | 6     | b_7 + b_11 |
    # |  3 vs 1   | 9     | b_7 + b_12 |
    # |  3 vs 1   | 12    | b_7 + b_13 |
    
    d_post_delta <- data.table(
      # 2 vs 1 at timepoint 1
      delta_2_1_1 = as.numeric(m_post[, 6]),
      # 2 vs 1 at timepoint 2
      delta_2_1_2 = as.numeric(m_post[, 6] + m_post[, 8]),
      # 2 vs 1 at timepoint 3 etc
      delta_2_1_3 = as.numeric(m_post[, 6] + m_post[, 9]),
      delta_2_1_4 = as.numeric(m_post[, 6] + m_post[, 10]),
      
      delta_3_1_1 = as.numeric(m_post[, 7]),
      delta_3_1_2 = as.numeric(m_post[, 7] + m_post[, 11]),
      delta_3_1_3 = as.numeric(m_post[, 7] + m_post[, 12]),
      delta_3_1_4 = as.numeric(m_post[, 7] + m_post[, 13])
    )
    
    # bind the intercept with betas
    d_post <- data.table(f_1$draws(variables = c("b_0"), format = "matrix"), m_post)
    names(d_post) <- paste0("b_", 0:13)
    d_post <- cbind(d_post, d_post_delta)
   
    # d_fig_1 <- melt(d_post, measure.vars = names(d_post))
    # 
    # ggplot(d_fig_1, aes(x = value)) +
    #   geom_density() +
    #   ggh4x::facet_wrap2(~variable, ncol = 2, scales = "free_x")
    
    if(return_posterior){
      
      
      d_tmp <- data.table(f_1$draws(variables = c("b_0", "b", "sigma", "rho"), format = "matrix"))
      d_post_all <- rbind(
        d_post_all,
        cbind(ic = l_spec$ic, d_tmp)
      )
      
      f_2 <- nlme::gls(y ~ y0_std + log(age0) + t_obs*trt, d_mod,
                         correlation = nlme::corAR1(form = ~ 1 | id))
      s <- summary(f_2)
      
      d_tmp <- data.table(par = rownames(s$tTable), value = s$coef, confint(f_2))
      d_tmp[, par_b := par_map[par]]
      
      # extract residual sd and autocorrel
      aux <- exp(as.numeric(f_2$modelStruct$corStruct))
      d_tmp <- rbind(
        d_tmp, 
        data.table(par = "se", value = s$sigma),
        data.table(par = "rho", value = (aux - 1) / (aux + 1)),
        fill = T
      )
      d_tmp[par == "se", par_b := "sigma"]
      d_tmp[par == "rho", par_b := "rho"]
      
      d_freq <- rbind(
        d_freq,
        cbind(ic = l_spec$ic, d_tmp)
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
        d_post_long[par %in% g_par_delta_pri, .(
          ic = l_spec$ic,
          rule = factor("ni", levels = g_rule_type),
          p = mean(value > l_spec$delta$ni),
          dec = as.integer(mean(value > l_spec$delta$ni) > l_spec$thresh$ni)
        ), keyby = par],
        d_post_long[par %in% g_par_delta_pri, .(
          ic = l_spec$ic,
          rule = factor("inf", levels = g_rule_type),
          p = mean(value < l_spec$delta$inf),
          dec = as.integer(mean(value < l_spec$delta$inf) > l_spec$thresh$inf)
        ), keyby = par]
      ),
      
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec  
      )
    ]
    
    log_info("Trial ", ix, " compared to thresholds ", l_spec$ic)
    
    # For trial stopping we only consider the comparisons to the
    # soc. In order to stop the study both need to have been resolved (either 
    # NI or inferiority have been concluded).
    d_stop <- d_pr_dec[
      ic <= l_spec$ic & par %in% c(g_par_delta_pri), 
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]
    
    
    # Update allocation probabilities if decisions have been reached.
    # Subsequent enrolments get redirected to the remaining arms.
    if(all(d_stop$resolved)){
      log_info("Stop trial as all questions addressed ", ix)
      stop_enrol <- T    
    } else if(any(d_stop$resolved)){
      
      # 2 vs 1 at timepoint 4    
      if(d_stop[par == "delta_2_1_4", resolved]) {
        l_spec$p_trt_alloc[2] <- 0  
        l_spec$p_trt_alloc <- l_spec$p_trt_alloc / sum(l_spec$p_trt_alloc)
      } 
      
      # 3 vs 1 at timepoint 4
      if(d_stop[par == "delta_3_1_4", resolved]) {
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
    # data collected in the trial
    
    # ic - enrolment cohort
    # ia - analysis in which pt is first analysed
    # id - pt id
    # trt - trt alloc
    # t0 - enrol day
    # tfu - fu day
    # fu - stage (0 = pre, 1 = post)
    # y - ppfev1
    # y_mis - missingness indicator
    # t_anlys - day on pt first enters their first analysis (early enrolments
    # will enter multiple analyses)
    d_all = d_all,
    
    d_post_smry_1 = d_post_smry_1,
    
    d_pr_dec = d_pr_dec,
    
    stop_at = stop_at
  )
  
  if(return_posterior){
    l_ret$d_post_all <- copy(d_post_all)
    l_ret$d_freq <- copy(d_freq)
  }
  # 
  
  return(l_ret)
}






run_sim03 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 3
  }
  
  l_spec <- list()
  
  l_spec$desc <- g_cfgsc$desc
  l_spec$nsim <- g_cfgsc$nsim
  l_spec$nex <- g_cfgsc$nex
  
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  
  l_spec$pt_per_day <-  g_cfgsc$pt_per_day
  l_spec$ramp_up_days <-  g_cfgsc$ramp_up_days
  
  l_spec$age_mu  <- g_cfgsc$age_mu
  l_spec$age_sig <- g_cfgsc$age_sig
  l_spec$age_lwr <- g_cfgsc$age_lwr
  l_spec$age_upr <- g_cfgsc$age_upr
  
  # missingness
  l_spec$pr_ymis <- g_cfgsc$pr_ymis
  
  # trt alloc
  l_spec$p_trt_alloc <- unlist(g_cfgsc$trt)/length(unlist(g_cfgsc$trt))
  
  l_spec$b_0 <- g_cfgsc$b_0
  l_spec$b_time <- unlist(g_cfgsc$b_time)
  l_spec$b_age <- g_cfgsc$b_age
  l_spec$b_trt <- do.call(rbind, g_cfgsc$b_trt)
  
  # subject levell heterogeneity - allows us to represent true (latent) mu0
  # and then what we actually observe (y0) 
  # data is generated based on mu0 and modelled adjusted for observed y0 (noisey
  # version of mu0)
  l_spec$sigma_i <- g_cfgsc$sigma_i
  
  l_spec$sigma <- g_cfgsc$sigma 
  l_spec$rho <- g_cfgsc$rho
  
  l_spec$t_sprty_obs <- unlist(g_cfgsc$t_sprty_obs)
  
  # # priors on log-odds/log-or
  # l_spec$prior <- list()
  # # location, scale
  # l_spec$prior$pri_b_0 <- unlist(g_cfgsc$pri_b_0)
  # l_spec$prior$pri_b_pre <- unlist(g_cfgsc$pri_b_pre)
  # l_spec$prior$pri_b_trt <- unlist(g_cfgsc$pri_b_trt)
  # l_spec$prior$pri_b_se <- g_cfgsc$pri_b_se
  # 
  # l_spec$prior$pri_b_0 <- unlist(g_cfgsc$pri_b_0)
  # l_spec$prior$pri_b_mu <- unlist(g_cfgsc$pri_b_mu)
  # l_spec$prior$pri_b_sd <- unlist(g_cfgsc$pri_b_sd)
  # l_spec$prior$pri_b_se <- g_cfgsc$pri_b_se
 
  # decision hurdle
  l_spec$delta <- list()
  l_spec$delta$ni <- g_cfgsc$dec_delta_ni
  l_spec$delta$inf <- g_cfgsc$dec_delta_inf
  
  # evidentiary requirement
  l_spec$thresh <- list()
  l_spec$thresh$ni <- unlist(g_cfgsc$dec_thresh_ni)
  l_spec$thresh$inf <- unlist(g_cfgsc$dec_thresh_inf)
  
  if(l_spec$nex > 0){
    log_info("Creating ", l_spec$nex, " example trials with full posterior")
    l_spec$ex_trial_ix <- sort(sample(1:g_cfgsc$nsim, size = l_spec$nex, replace = F))
  }
  return_posterior <- F
  str(l_spec)
  e = NULL
  ix <- 1
  
  ## LOOP -------
  
  
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
  
  d_freq <- data.table(do.call(rbind, lapply(1:length(r), function(i){
    # if the sim contains full posterior (for example trial) then return
    if(!is.null(r[[i]]$d_freq)){
      cbind(sim = i, r[[i]]$d_freq)
    }
    
  } )))
  
  l <- list(
    cfg = l_spec,
    d_pr_dec = d_pr_dec, 
    d_post_smry_1 = d_post_smry_1,
    d_all = d_all,
    d_post_all = d_post_all,
    d_freq = d_freq
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim03/sim03-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim03 <- function(){
  log_info("run_none_sim03: Nothing doing here bud.")
}

main_sim03 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim03()


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
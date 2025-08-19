library(data.table)






get_sim01_trial_data <- function(
    l_spec
){
  
  if(is.null(l_spec$ia)){
    ia <- 1
  } else {
    ia <- l_spec$ia
  }
  if(is.null(l_spec$is)){
    is <- 1
  } else {
    is <- l_spec$is
  }
  if(is.null(l_spec$ie)){
    ie <- 1
  } else {
    ie <- l_spec$ie
  }
  if(is.null(l_spec$t0)){
    t0 <- 1
  } else {
    t0 <- l_spec$t0
  }
  
  
  d <- data.table(
    ia = ia,
    id = is:ie,
    t0 = t0,
    lung = rbinom(l_spec$N[ia], 1, prob = l_spec$p_lung_alloc) + 1
  )
  # lung 40%  (< 70% pFEV1) = 1, (>= 70% pFEV1) = 2
  
  # colonisation (no pseudo) = 1, (pseudo) = 2 
  d[lung == 1, pseudo := rbinom(.N, 1, prob = l_spec$p_pseudo_alloc_lt70) + 1]
  d[lung == 2, pseudo := rbinom(.N, 1, prob = l_spec$p_pseudo_alloc_ge70) + 1]
  
  # lung (< 70% pFEV1), no pseudo
  trt_opts <- 1:length(l_spec$p_trt_alloc)
  
  d[lung == 1 & pseudo == 1, trt := sample(trt_opts, size = .N, replace = T, prob = l_spec$p_trt_alloc)]
  # lung (< 70% pFEV1), pseudo
  d[lung == 1 & pseudo == 2, trt := sample(trt_opts, size = .N, replace = T, prob = l_spec$p_trt_alloc)]
  # lung (>= 70% pFEV1), no pseudo
  d[lung == 2 & pseudo == 1, trt := sample(trt_opts, size = .N, replace = T, prob = l_spec$p_trt_alloc)]
  # lung (>= 70% pFEV1), pseudo
  d[lung == 2 & pseudo == 2, trt := sample(trt_opts, size = .N, replace = T, prob = l_spec$p_trt_alloc)]
  
  # just for sanity checking
  
  # d_tbl <- d[, .N/nrow(d), keyby = .(lung, pseudo)]
  # d_tbl[, lung := factor(lung, labels = c("< 70% pFEV1", ">= 70% pFEV1"), levels = 1:2)]
  # d_tbl[, pseudo := factor(pseudo, labels = c("no pseudo", "pseudo"), levels = 1:2)]
  # pander::pandoc.table(d_tbl)
  
  # d_tbl <- d[, .N, keyby = .(lung, pseudo, trt)]
  # pander::pandoc.table(d_tbl)
  
  # Given the simplicity of the model we can specify the linear predictor 
  # directly in terms of risk increments.
  d[, p := l_spec$bmu + l_spec$blung[lung] + l_spec$bpseudo[pseudo] + l_spec$btrt[trt]]
  # y = 1 indicates treatment failure
  d[, y := rbinom(.N, 1, p)]
  
  d
}



get_sim01_stan_data <- function(d_all){
  
  # convert from binary representation to binomial (successes/trials)
  d_mod <- d_all[, .(y = sum(y), n = .N, eta = unique(p)), 
                 keyby = .(lung, pseudo, trt)]
  
  d_mod[, p_obs := y / n]
  
  ld <- list(
    # full dataset
    N = nrow(d_mod), 
    y = d_mod[, y], 
    n = d_mod[, n], 
    trt = d_mod[, trt], 
    lung = d_mod[, lung], 
    pseudo = d_mod[, pseudo],
    prior_only = 0
  )
  
  list(
    d_mod = d_mod,
    ld = ld
  )
  
}
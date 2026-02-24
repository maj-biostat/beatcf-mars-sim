

library(data.table)
library(pbapply)
library(survival)
library(alda)
library(cmdstanr)
library(ggplot2)
library(truncnorm)
library(lme4)
library(splines)
library(INLA)

logit_inv <- function(x) 1 / (1 + exp(-x))


baseline_ppfev <- function(age) {
  100 - 45 * (log(age) - log(10)) / (log(60) - log(10))
}


hazard_hex <- function(t_health, ppfev) {
  linpred <-
    -5 +                  # baseline (tuned to ~2 exacerbations/year)
    0.001 * t_health +       # increasing hazard with duration healthy
    -0.001 * (ppfev - 80)    # worse lung function -> higher hazard
  
  logit_inv(linpred)
}

hazard_rec <- function(t_exac, ppfev, trt) {
  
  trt_effect <- c(
    soc     = 0,
    delayed = -0.08,  # slower recovery
    deferred = -0.08
  )[trt]
  
  linpred <-
    -3.5 +                    # baseline
    0.2 * t_exac +           # strong time dependence
    -0.04 * (ppfev - 80) +    # sicker recover slower
    trt_effect
  
  logit_inv(linpred)
}



n_patients <- 300
n_days <- 365

patients <- data.table(
  id = 1:n_patients,
  age = rtruncnorm(
    n_patients,
    a = log(10),
    b = log(60),
    mean = log(30),
    sd = 0.5
  )
)

patients[, age := exp(age)]
patients[, ppfev_base := baseline_ppfev(age)]


sample_trt <- function() {
  sample(c("soc", "delayed", "deferred"), 1)
}



sim_list <- vector("list", n_patients)

i <- 1
for (i in seq_len(n_patients)) {
  
  age <- patients[id == i, age]
  ppfev_base <- patients[id == i, ppfev_base]
  
  dt <- data.table(
    id = i,
    day = 1:n_days,
    state = NA_integer_,
    ppfev = NA_real_,
    t_state = NA_integer_,
    trt = NA_character_
  )
  
  # initial conditions
  state <- 0
  t_state <- 1
  trt <- NA_character_
  
  ppfev <- ppfev_base
  
  ppfev_pre_exac <- ppfev_base
  
  d <- 1
  for (d in 1:n_days) {
    
    # secular decline ~1% per year
    ppfev <- ppfev * exp(-0.01 / 365)
    
    if (state == 0) {
      
      haz <- hazard_hex(t_state, ppfev)
      event <- rbinom(1, 1, haz)
      
      if (event == 1) {
        state <- 1
        t_state <- 1
        trt <- sample_trt()
        
        ppfev_pre_exac <- ppfev
        ppfev <- ppfev * (1 - 0.03 + rnorm(1, 0, 0.05))
      } else {
        t_state <- t_state + 1
      }
      
    } else {
      
      haz <- hazard_rec(t_state, ppfev, trt)
      event <- rbinom(1, 1, haz)
      
      # recovery dynamics
      target <- ppfev_pre_exac
      ppfev <- ppfev + 0.15 * (target - ppfev)
      
      if (event == 1) {
        state <- 0
        t_state <- 1
        trt <- NA_character_
      } else {
        t_state <- t_state + 1
      }
    }
    
    # noise
    ppfev <- ppfev + rnorm(1, 0, 0.3)
    
    set(dt, i = which(dt$day == d), j = "state", state)
    set(dt, i = which(dt$day == d), j = "ppfev", ppfev)
    set(dt, i = which(dt$day == d), j = "t_state", t_state)
    set(dt, i = which(dt$day == d), j = "trt", trt)
    
    # dt[day == d]$state <- state
    # dt[day == d]$ppfev <- ppfev
    # dt[day == d]$t_state <- t_state
    # dt[day == d]$trt <- trt
  }
  
  sim_list[[i]] <- dt
}


sim_dt <- rbindlist(sim_list)
setorder(sim_dt, id, day)

sim_dt[id == 1]

# run-length–encoding to move to start stop format:
d_w <- sim_dt[, run := rleid(state)][, .(
  start = min(day),
  stop  = max(day),
  state = state[1L]
),
by = .(id, run)]
sim_dt[, run := NULL]
d_w[, dur := stop - start + 1]
d_w[id == 1]

# how many have pex in year? 
d_w[, .N, keyby = .(id, state)][state == 1, .N]
# of those that have pex, how many do they have in the year?
d_w[, .N, keyby = .(id, state)][state == 1, .(mean(N), min(N), max(N))]
# of those who have pex, how long are they in days?
d_w[state == 1, .(mean(dur), min(dur), max(dur))]


# ggplot(sim_dt[id %in% sample(1:n_patients, size = 10, replace = F)], 
#        aes(x = day, y = ppfev, group = id)) +
#   geom_line(lwd = 0.3) +
#   geom_point(aes(col = factor(state)), size = 0.2) +
#   theme(legend.position = "bottom")
#



# You have two transitions:
#   
#   0 → 1 : Healthy → Exacerbation
# 
# 1 → 0 : Exacerbation → Recovery
# 
# This is a two-state recurrent semi-Markov model because transition probabilities depend on:
#   
#   time in current state (t_state)
# 
# time-varying ppFEV
# 
# treatment (for recovery)
# 
# 
# You need to construct a person-day risk set with a transition indicator.
# 
# From your simulated dataset sim_dt, construct:
#   
#   event01 = 1 if transition 0→1 at day d
# 
# event10 = 1 if transition 1→0 at day d

dt <- copy(sim_dt)


dt[, prev_state := shift(state), by = id]
dt <- dt[!is.na(prev_state)]

dt[, trans := paste0(prev_state, state)]
dt[, trans := factor(trans, levels = c("00","01","10","11"))]

dt[, prev_state_f := factor(prev_state)]

library(nnet)

fit_multi <- multinom(
  trans ~ 
    prev_state_f +
    prev_state_f:t_state +
    prev_state_f:I(t_state^2) +
    prev_state_f:ppfev +
    prev_state_f:factor(trt),
  data = dt,
  trace = FALSE
)


# only interested in those that have an exacerbation
dt <- dt[id %in% dt[state == 1, unique(id)]]

dt[, prev_state := shift(state), by = id]

# safe as first day is fixed state
dt[, event01 := as.integer(prev_state == 0 & state == 1)]
dt[, event10 := as.integer(prev_state == 1 & state == 0)]

# Remove first day per patient:
dt <- dt[!is.na(prev_state)]

View(dt[prev_state == 0])


fit01_re <- glmer(
  event01 ~ factor(t_state) + ppfev + (1 | id),
  data = dt[prev_state == 0],
  family = binomial()
)

fit01 <- glm(event01 ~ t_state + I(t_state^2) + ppfev,
    family = binomial(),
    data = dt[prev_state == 0])

# fit01_re <- inla(
#   event01 ~ ns(t_state,4) + ppfev + f(id, model = "iid"), 
#   data = dt[prev_state == 0],
#   family = binomial()
# )
# fit01_re <- glm(
#   event01 ~ 
#     ns(t_state, df = 4) + 
#     ppfev,
#   data = dt[prev_state == 0],
#   family = binomial()
# )

fit10_re <- glmer(
  event10 ~ factor(t_state) + ppfev + factor(trt) + (1 | id),
  data = dt[prev_state == 1],
  family = binomial()
)









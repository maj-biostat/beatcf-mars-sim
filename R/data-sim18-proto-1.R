# uses matrices of probabilities to control the state transitions.
# this means we know the underlying truth of the target we are trying to 
# estimate.

library(data.table)
library(ggplot2)

##------------------------------------------------------------
## Utilities
##------------------------------------------------------------

inv_logit <- function(x)
  1 / (1 + exp(-x))

## Draw from proportional odds model
rord_pom <- function(lp, alpha) {
  
  p1 <- inv_logit(alpha[1] - lp)
  p2 <- inv_logit(alpha[2] - lp) - p1
  p3 <- 1 - p1 - p2
  
  sample.int(
    3L,
    size = 1L,
    prob = c(p1, p2, p3)
  ) 
}

##------------------------------------------------------------
## Transition probability model
##
## This is the ONLY part that should normally be edited.
##
## Returns probabilities for
##
## 0 = None
## 1 = Mild
## 2 = Severe
##
## conditional on previous state, treatment and day.
##------------------------------------------------------------

transition_probs <- function(prev,
                             trt,
                             day)
{
  
  ## baseline recovery probability
  ##
  ## increases rapidly over first week then plateaus
  
  r <- plogis(-2 + 0.35 * day)
  
  ## treatment effects
  ##
  ## positive = slower recovery
  
  if (trt == "soc")
    tx <- 1.00
  
  if (trt == "discon")
    tx <- 0.92
  
  if (trt == "def")
    tx <- 0.80
  
  r <- r * tx
  
  ##--------------------------------------------------------
  ## Previous NONE
  ##--------------------------------------------------------
  
  if (prev == 1) {
    
    p <- c(
      0.985,
      0.013,
      0.002
    )
    
  }
  
  ##--------------------------------------------------------
  ## Previous MILD
  ##--------------------------------------------------------
  
  if (prev == 2) {
    
    p_none <- r
    
    p_severe <- (1 - r) * 0.15
    
    p_mild <- 1 - p_none - p_severe
    
    p <- c(
      p_none,
      p_mild,
      p_severe
    )
    
  }
  
  ##--------------------------------------------------------
  ## Previous SEVERE
  ##--------------------------------------------------------
  
  if (prev == 3) {
    
    p_none <- 0.15 * r
    
    p_mild <- 0.85 * r
    
    p_severe <- 1 - p_none - p_mild
    
    p <- c(
      p_none,
      p_mild,
      p_severe
    )
    
  }
  
  p / sum(p)
  
}

##------------------------------------------------------------
## Convert desired probabilities to proportional odds LP
##
## alpha is fixed.
##
## The returned LP produces approximately the supplied
## probabilities under the proportional odds model.
##------------------------------------------------------------

approx_lp <- function(prob,
                      alpha)
{
  
  F1 <- prob[1]
  
  F2 <- prob[1] + prob[2]
  
  eta1 <- alpha[1] - qlogis(F1)
  
  eta2 <- alpha[2] - qlogis(F2)
  
  mean(c(eta1, eta2))
  
}

##------------------------------------------------------------
## Simulator
##------------------------------------------------------------
# days = 28
# visit_days = c(0:14, 21, 28)
# alpha = c(-1, 1)
sim_cf_trial <- function(
    n,
    days = 28,
    visit_days = c(0:14, 21, 28),
    alpha = c(-1, 1)
    )
{
  
  trt <- sample(
    c("soc", "discon", "def"),
    size = n,
    replace = TRUE
  )
  
  dt <- CJ(
    id = seq_len(n),
    day = 0:days
  )
  
  dt[, trt := trt[id]]
  
  setorder(dt, id, day)
  
  dt[, state := NA_integer_]
  
  # just assume 50/50 mild/sev
  dt[day == 0, state := sample(2:3, .N, replace = T)]
  
  ##--------------------------------------------------------
  ## Generate latent daily states
  ##--------------------------------------------------------
  
  for (ii in seq_len(n)) {
    
    rows <- which(dt$id == ii)
    
    for (j in 2:length(rows)) {
      
      prev <- dt$state[rows[j - 1]]
      
      day <- dt$day[rows[j]]
      
      trt <- dt$trt[rows[j]]
      
      ## Clinically interpretable transition probabilities
      
      p_target <- transition_probs(
        prev = prev,
        trt = trt,
        day = day
      )
      
      ## Corresponding proportional odds LP
      
      lp <- approx_lp(
        prob = p_target,
        alpha = alpha
      )
      
      ## Draw next state from proportional odds model
      
      dt$state[rows[j]] <-
        rord_pom(
          lp = lp,
          alpha = alpha
        )
      
    }
    
  }
  
  ##--------------------------------------------------------
  ## Observation schedule
  ##--------------------------------------------------------
  
  obs <- dt[
    day %in% visit_days
  ]
  
  ##--------------------------------------------------------
  ## Expected symptom-free days
  ##
  ## (using latent daily process)
  ##--------------------------------------------------------
  
  symptom_free <- dt[
    ,
    .(
      symptom_free_days = sum(state == 1)
    ),
    by = .(id, trt)
  ]
  
  ##--------------------------------------------------------
  ## State occupancy
  ##--------------------------------------------------------
  
  occupancy <- dt[
    ,
    .N,
    by = .(
      trt,
      day,
      state
    )
  ]
  
  occupancy[
    ,
    prob := N / sum(N),
    by = .(trt, day)
  ]
  
  ##--------------------------------------------------------
  ## Recrudescence
  ##
  ## Simple definition:
  ##
  ## once recovered,
  ## subsequently develops symptoms.
  ##
  ##--------------------------------------------------------
  
  recrudescence <- dt[
    ,
    {
      
      recovered <- which(state == 1)
      
      recur <- FALSE
      
      if (length(recovered) > 0) {
        
        first <- recovered[1]
        
        recur <- any(
          state[(first + 1):.N] > 1
        )
        
      }
      
      .(
        recrudescence = recur
      )
      
    },
    by = .(
      id,
      trt
    )
  ]
  
  list(
    
    daily = dt,
    
    observed = obs,
    
    symptom_free = symptom_free,
    
    occupancy = occupancy,
    
    recrudescence = recrudescence
    
  )
  
}

##------------------------------------------------------------
## Example
##------------------------------------------------------------

set.seed(123)

sim <- sim_cf_trial(
  n = 600
)

## Mean symptom-free days

sim$symptom_free[
  ,
  .(
    mean_days = mean(symptom_free_days),
    sd_days = sd(symptom_free_days)
  ),
  by = trt
]

## State occupancy

sim$occupancy[order(state, day)][trt == "soc"]

## Recrudescence rates

sim$recrudescence[
  ,
  .(
    recrudescence = mean(recrudescence)
  ),
  by = trt
]

## Observed data (what would actually be collected)

head(sim$observed)










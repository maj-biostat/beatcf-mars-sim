# Simulation Study for Operating Characteristics
# Evaluate power to detect non-inferiority under different scenarios

library(data.table)
library(survival)
library(coxme)
library(lme4)
library(lmerTest)
library(parallel)
library(emmeans)

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

# Study design
n_patients <- 600
n_days <- 365
n_simulations <- 10  # Start with 100, increase for final analysis

# True treatment effects to simulate under
scenarios <- list(
  # Scenario 1: Null (no difference)
  list(name = "Null", 
       hr_delayed = 1.0, hr_deferred = 1.0,
       ppfev_diff_delayed = 0, ppfev_diff_deferred = 0),
  
  # Scenario 2: Small inferior effect (at NI margin)
  list(name = "At_NI_margin",
       hr_delayed = 0.7, hr_deferred = 0.7,
       ppfev_diff_delayed = -3, ppfev_diff_deferred = -3),
  
  # Scenario 3: Moderate inferior effect (beyond NI margin)
  list(name = "Beyond_NI_margin",
       hr_delayed = 0.6, hr_deferred = 0.6,
       ppfev_diff_delayed = -5, ppfev_diff_deferred = -5),
  
  # Scenario 4: Differential effects (delayed better than deferred)
  list(name = "Differential",
       hr_delayed = 0.85, hr_deferred = 0.7,
       ppfev_diff_delayed = -1, ppfev_diff_deferred = -3)
)

# Non-inferiority margins
ni_margin_hr <- 0.7
ni_margin_ppfev <- -3

cat("SIMULATION STUDY: OPERATING CHARACTERISTICS\n")


cat("Study parameters:\n")
cat("  Sample size:", n_patients, "patients\n")
cat("  Follow-up:", n_days, "days\n")
cat("  Simulations per scenario:", n_simulations, "\n")
cat("  NI margin (HR):", ni_margin_hr, "\n")
cat("  NI margin (ppFEV):", ni_margin_ppfev, "%\n\n")

# ============================================================================
# DATA GENERATION FUNCTION
# ============================================================================

simulate_trial_dataset <- function(seed, hr_delayed, hr_deferred, 
                                   ppfev_diff_delayed, ppfev_diff_deferred) {
  set.seed(seed)
  
  # Generate patients
  ages <- rlnorm(n_patients, meanlog = log(35), sdlog = 0.4)
  ages <- pmin(pmax(ages, 10), 60)
  initial_ppfevs <- 100 - (ages - 10) * (45 / 50)
  
  daily_decline_healthy <- 1 / 365
  exacerbation_drop_mean <- 3
  exacerbation_drop_sd <- 0.5
  recovery_half_life <- 5
  
  # Functions with treatment effect parameters
  calc_exacerbation_hazard <- function(duration_healthy, ppfev) {
    ppfev_scaled <- (ppfev - 80) / 20
    log_hazard <- -6.5 + 0.3 * log(duration_healthy + 1) - 0.5 * ppfev_scaled
    plogis(log_hazard)
  }
  
  calc_recovery_hazard <- function(duration_exacerb, ppfev, treatment, 
                                   hr_delayed, hr_deferred) {
    ppfev_scaled <- (ppfev - 80) / 20
    
    # Convert HR to log-HR for additive effect
    treatment_effect <- 0
    if (treatment == "delayed") {
      treatment_effect <- log(hr_delayed)
    } else if (treatment == "deferred") {
      treatment_effect <- log(hr_deferred)
    }
    
    log_hazard <- -1.5 + 0.4 * duration_exacerb + 0.4 * ppfev_scaled + treatment_effect
    plogis(log_hazard)
  }
  
  # Simulate patients
  all_exacerbations <- list()
  all_daily <- list()
  
  for (i in 1:n_patients) {
    state <- "healthy"
    duration_in_state <- 0
    ppfev <- initial_ppfevs[i]
    ppfev_before_exacerbation <- initial_ppfevs[i]
    exacerbation_drop <- 0
    current_treatment <- NA
    exacerbation_count <- 0
    
    for (day in 1:n_days) {
      duration_in_state <- duration_in_state + 1
      
      if (state == "healthy") {
        ppfev <- ppfev - daily_decline_healthy
        
        p_exacerbation <- calc_exacerbation_hazard(duration_in_state, ppfev)
        if (runif(1) < p_exacerbation) {
          state <- "exacerbation"
          duration_in_state <- 0
          ppfev_before_exacerbation <- ppfev
          exacerbation_count <- exacerbation_count + 1
          
          current_treatment <- sample(c("standard", "delayed", "deferred"), 1)
          
          exacerbation_drop <- rnorm(1, exacerbation_drop_mean, exacerbation_drop_sd)
          exacerbation_drop <- max(exacerbation_drop, 1)
          
          # Add treatment-specific ppFEV effect at onset
          if (current_treatment == "delayed") {
            exacerbation_drop <- exacerbation_drop - ppfev_diff_delayed
          } else if (current_treatment == "deferred") {
            exacerbation_drop <- exacerbation_drop - ppfev_diff_deferred
          }
          
          ppfev <- ppfev - exacerbation_drop
          episode_start_day <- day
          episode_start_ppfev <- ppfev
        }
        
      } else if (state == "exacerbation") {
        recovery_rate <- log(2) / recovery_half_life
        daily_recovery <- exacerbation_drop * recovery_rate * 
          exp(-recovery_rate * duration_in_state)
        ppfev <- min(ppfev + daily_recovery, ppfev_before_exacerbation)
        
        p_recovery <- calc_recovery_hazard(duration_in_state, ppfev, 
                                           current_treatment,
                                           hr_delayed, hr_deferred)
        if (runif(1) < p_recovery) {
          episode_duration <- duration_in_state + 1
          
          all_exacerbations[[length(all_exacerbations) + 1]] <- data.table(
            patient_id = i,
            episode_number = exacerbation_count,
            start_day = episode_start_day,
            end_day = day,
            duration = episode_duration,
            treatment = current_treatment,
            ppfev_at_onset = episode_start_ppfev,
            ppfev_at_recovery = ppfev,
            censored = FALSE
          )
          
          state <- "healthy"
          duration_in_state <- 0
          exacerbation_drop <- 0
          current_treatment <- NA
        }
      }
      
      all_daily[[length(all_daily) + 1]] <- data.table(
        patient_id = i,
        day = day,
        state = state,
        treatment = ifelse(is.na(current_treatment), "none", current_treatment),
        ppfev = ppfev,
        exacerbation_number = exacerbation_count
      )
    }
  }
  
  df_exacerb <- rbindlist(all_exacerbations)
  df_exacerb[, treatment := factor(treatment, levels = c("standard", "delayed", "deferred"))]
  df_daily <- rbindlist(all_daily)
  
  list(exacerbations = df_exacerb, daily = df_daily)
}

# ============================================================================
# ANALYSIS FUNCTION
# ============================================================================

analyze_dataset <- function(data) {
  # Extract exacerbation data
  exacerb <- data$exacerbations
  
  if (nrow(exacerb) < 50) {
    # Not enough data
    return(list(
      recovery_delayed_ni = NA,
      recovery_deferred_ni = NA,
      ppfev_delayed_ni = NA,
      ppfev_deferred_ni = NA,
      n_exacerb = nrow(exacerb)
    ))
  }
  
  # Analysis 1: Cox model for recovery
  tryCatch({
    cox_model <- coxme(Surv(duration, rep(1, len = nrow(exacerb))) ~ treatment + (1 | patient_id),
                       data = exacerb)
    
    cox_coef <- summary(cox_model)$coefficients
    
    # Extract treatment effects
    beta_delayed <- cox_coef["treatmentdelayed", "coef"]
    se_delayed <- cox_coef["treatmentdelayed", "se(coef)"]
    beta_deferred <- cox_coef["treatmentdeferred", "coef"]
    se_deferred <- cox_coef["treatmentdeferred", "se(coef)"]
    
    # NI test: lower CI of log(HR) > log(0.7)
    log_ni_margin <- log(ni_margin_hr)
    ci_delayed_lower <- beta_delayed - 1.96 * se_delayed
    ci_deferred_lower <- beta_deferred - 1.96 * se_deferred
    
    recovery_delayed_ni <- ci_delayed_lower > log_ni_margin
    recovery_deferred_ni <- ci_deferred_lower > log_ni_margin
    
  }, error = function(e) {
    recovery_delayed_ni <- NA
    recovery_deferred_ni <- NA
  })
  
  # Analysis 2: ppFEV at day 14
  # Prepare longitudinal data
  daily <- data$daily
  
  exacerb_long_list <- list()
  i <- 1
  for (i in 1:nrow(exacerb)) {
    pid <- exacerb$patient_id[i]
    start_day <- exacerb$start_day[i]
    end_day <- exacerb$end_day[i]
    episode_num <- exacerb$episode_number[i]
    treatment <- exacerb$treatment[i]
    
    max_day <- min(start_day + 14, end_day)
    
    episode_days <- daily[patient_id == pid & day >= start_day & day <= max_day]
    episode_days[, time_from_onset := day - start_day]
    episode_days[, episode_id := paste0(pid, "_", episode_num)]
    episode_days[, episode_treatment := treatment]
    
    exacerb_long_list[[i]] <- episode_days
  }
  
  exacerb_long <- rbindlist(exacerb_long_list, fill = TRUE)
  
  tryCatch({
    lmm_model <- lmer(ppfev ~ episode_treatment * time_from_onset + 
                        (1 | patient_id) + (1 | episode_id),
                      data = exacerb_long)
    
    # Get contrasts at day 14
    emm_day14 <- suppressMessages(emmeans(lmm_model, ~ episode_treatment, 
                         at = list(time_from_onset = 14)))
    contrasts_day14 <- contrast(emm_day14, method = "trt.vs.ctrl", ref = "standard")
    contrasts_summary <- summary(contrasts_day14, infer = TRUE)
    
    # NI test: lower CI > -3%
    ppfev_delayed_ni <- contrasts_summary$asymp.LCL[
      contrasts_summary$contrast == "delayed - standard"] > ni_margin_ppfev
    ppfev_deferred_ni <- contrasts_summary$asymp.LCL[
      contrasts_summary$contrast == "deferred - standard"] > ni_margin_ppfev
    
  }, error = function(e) {
    ppfev_delayed_ni <- NA
    ppfev_deferred_ni <- NA
  })
  
  list(
    recovery_delayed_ni = recovery_delayed_ni,
    recovery_deferred_ni = recovery_deferred_ni,
    ppfev_delayed_ni = ppfev_delayed_ni,
    ppfev_deferred_ni = ppfev_deferred_ni,
    n_exacerb = nrow(exacerb)
  )
}

# ============================================================================
# RUN SIMULATIONS
# ============================================================================

results_all <- list()

scenario <- scenarios[[1]]
sim_id <- 1

# for (scenario in scenarios) {
  
  scenario <- scenarios[[1]]
  
  cat("\n")
  cat("SCENARIO:", scenario$name, "\n")
  cat("  HR delayed:", scenario$hr_delayed, "\n")
  cat("  HR deferred:", scenario$hr_deferred, "\n")
  cat("  ppFEV diff delayed:", scenario$ppfev_diff_delayed, "%\n")
  cat("  ppFEV diff deferred:", scenario$ppfev_diff_deferred, "%\n\n")
  
  cat("Running", n_simulations, "simulations...\n")
  
  # Run simulations in parallel
  n_cores <- min(detectCores() - 1, 4)
  
  sim_results <- mclapply(1:n_simulations, function(sim_id) {
    if (sim_id %% 10 == 0) {
      cat("  Simulation", sim_id, "\n")
    }
    
    # tictoc::tic()
    # Generate data
    data <- simulate_trial_dataset(
      seed = sim_id + 1000,
      hr_delayed = scenario$hr_delayed,
      hr_deferred = scenario$hr_deferred,
      ppfev_diff_delayed = scenario$ppfev_diff_delayed,
      ppfev_diff_deferred = scenario$ppfev_diff_deferred
    )
    
    # tictoc::toc()
    # Analyze
    analyze_dataset(data)
    
    
  }, mc.cores = n_cores)
  
  # Combine results
  results_dt <- rbindlist(sim_results)
  
  # Calculate operating characteristics
  cat("\nOperating Characteristics:\n")
  cat("-"*70, "\n\n")
  
  cat("RECOVERY (Time to recovery - Cox model):\n")
  cat(sprintf("  Delayed - Non-inferiority power: %.1f%%\n",
              100 * mean(results_dt$recovery_delayed_ni, na.rm = TRUE)))
  cat(sprintf("  Deferred - Non-inferiority power: %.1f%%\n\n",
              100 * mean(results_dt$recovery_deferred_ni, na.rm = TRUE)))
  
  cat("ppFEV AT DAY 14 (Linear mixed model):\n")
  cat(sprintf("  Delayed - Non-inferiority power: %.1f%%\n",
              100 * mean(results_dt$ppfev_delayed_ni, na.rm = TRUE)))
  cat(sprintf("  Deferred - Non-inferiority power: %.1f%%\n\n",
              100 * mean(results_dt$ppfev_deferred_ni, na.rm = TRUE)))
  
  cat(sprintf("Mean exacerbations per dataset: %.1f\n", 
              mean(results_dt$n_exacerb, na.rm = TRUE)))
  
  # Store results
  results_all[[scenario$name]] <- results_dt
# }

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n")
cat("SUMMARY OF OPERATING CHARACTERISTICS\n")

# summary_table <- data.table(
#   Scenario = names(results_all),
#   Recovery_Delayed_Power = sapply(results_all, function(x) 
#     100 * mean(x$recovery_delayed_ni, na.rm = TRUE)),
#   Recovery_Deferred_Power = sapply(results_all, function(x) 
#     100 * mean(x$recovery_deferred_ni, na.rm = TRUE)),
#   ppFEV_Delayed_Power = sapply(results_all, function(x) 
#     100 * mean(x$ppfev_delayed_ni, na.rm = TRUE)),
#   ppFEV_Deferred_Power = sapply(results_all, function(x) 
#     100 * mean(x$ppfev_deferred_ni, na.rm = TRUE))
# )
# 
# print(summary_table)
# 
# # Save results
# save(results_all, summary_table, file = "simulation_results.RData")
# fwrite(summary_table, "simulation_summary.csv")
# 
# cat("\n\nResults saved to: simulation_results.RData and simulation_summary.csv\n")
# cat("\nSimulation study complete!\n")

# Simulation of Patient State and ppFEV Trajectories
# Using discrete-time continuation ratio models with data.table

library(data.table)
library(ggplot2)

set.seed(123)

# Parameters ----------------------------------------------------------------
n_patients <- 100
n_days <- 365

# Helper Functions ----------------------------------------------------------

# Age distribution (log-normal to get ages between ~10 and ~60)
generate_ages <- function(n) {
  ages <- rlnorm(n, meanlog = log(35), sdlog = 0.4)
  ages <- pmin(pmax(ages, 10), 60)  # Clip to [10, 60]
  return(ages)
}

# Assume y(x) = A exp (-k x^p) with p < 1 this gives steep early decline then
# flattens out (Weibull). Pick two points you want to joint together, e.g.
# f(10) == 100, f(60) = 55 and p <1 then solve for k and A
make_stretched_exp <- function(x1, y1, x2, y2, p = 0.6) {
  k <- -log(y2 / y1) / (x2^p - x1^p)
  A <- y1 * exp(k * x1^p)
  function(x) A * exp(-k * x^p)
}

f <- make_stretched_exp(10, 100, 60, 55, p = 0.3)

# Initial ppFEV based on age
# Just do linear decline - at age 10: ~100%, at age 60: ~55% 
initial_ppfev <- function(age) {
  f(age)
}

# Logistic function (inverse logit)
expit <- function(x) {
  1 / (1 + exp(-x))
}

# Calculate probability of transitioning from healthy to exacerbation
calc_exacerbation_hazard <- function(duration_healthy, ppfev) {
  ppfev_scaled <- (ppfev - 80) / 20  # Center and scale
  log_hazard <- -6.5 + 0.3 * log(duration_healthy + 1) - 0.5 * ppfev_scaled
  expit(log_hazard)
}

# Calculate probability of transitioning from exacerbation to healthy
calc_recovery_hazard <- function(duration_exacerb, ppfev, treatment) {
  ppfev_scaled <- (ppfev - 80) / 20
  
  # Treatment effect (delayed and deferred increase duration by 5%)
  treatment_effect <- ifelse(treatment %in% c("delayed", "deferred"), -0.3, 0)
  
  # Baseline hazard increases rapidly with duration
  log_hazard <- -2.5 + 0.3 * duration_exacerb + 0.4 * ppfev_scaled + treatment_effect
  expit(log_hazard)
}

# Initialize Patient Data ---------------------------------------------------
patients <- data.table(
  patient_id = 1:n_patients,
  age = generate_ages(n_patients)
)
patients[, initial_ppfev := initial_ppfev(age)]

# ggplot(patients, aes(x = age, y = initial_ppfev)) +
#   geom_point()


# Simulation Parameters
daily_decline_healthy <- 1 / 365  # 1% per year
daily_decline_healthy_sd <- 20 / 365

exacerbation_drop_mean <- 5
exacerbation_drop_sd <- 0.5
recovery_half_life <- 5  # days to recover 50%

# Simulate One Patient ------------------------------------------------------
simulate_patient <- function(patient_row) {
  pid <- patient_row$patient_id
  age <- patient_row$age
  baseline_ppfev <- patient_row$initial_ppfev
  
  # Initialize tracking variables
  state <- "healthy"
  duration_in_state <- 0
  ppfev <- baseline_ppfev
  ppfev_before_exacerbation <- baseline_ppfev
  exacerbation_drop <- 0
  current_treatment <- NA  # No treatment until exacerbation
  
  # Storage for daily observations
  daily_data <- vector("list", n_days)
  
  for (day in 1:n_days) {
    duration_in_state <- duration_in_state + 1
    
    # Update ppFEV based on current state
    if (state == "healthy") {
      # Gradual decline in healthy state
      ppfev <- ppfev - daily_decline_healthy + rnorm(1, 0, daily_decline_healthy_sd)
      
      # Check for transition to exacerbation
      p_exacerbation <- calc_exacerbation_hazard(duration_in_state, ppfev)
      if (runif(1) < p_exacerbation) {
        # Transition to exacerbation
        state <- "exacerbation"
        duration_in_state <- 0
        ppfev_before_exacerbation <- ppfev
        
        # RANDOMIZE TREATMENT AT EXACERBATION ONSET
        current_treatment <- sample(c("standard", "delayed", "deferred"), 1)
        
        # Apply exacerbation drop with noise
        exacerbation_drop <- rnorm(1, exacerbation_drop_mean, exacerbation_drop_sd)
        exacerbation_drop <- max(exacerbation_drop, 1)  # At least 1% drop
        ppfev <- ppfev - exacerbation_drop
      }
      
    } else if (state == "exacerbation") {
      # Gradual recovery during exacerbation (exponential)
      recovery_rate <- log(2) / recovery_half_life
      daily_recovery <- exacerbation_drop * recovery_rate * exp(-recovery_rate * duration_in_state)
      ppfev <- min(ppfev + daily_recovery, ppfev_before_exacerbation)
      
      # Check for transition to healthy
      p_recovery <- calc_recovery_hazard(duration_in_state, ppfev, current_treatment)
      if (runif(1) < p_recovery) {
        # Transition to healthy
        state <- "healthy"
        duration_in_state <- 0
        exacerbation_drop <- 0
        current_treatment <- NA  # Clear treatment after recovery
      }
    }
    
    # Store daily observation
    daily_data[[day]] <- data.table(
      patient_id = pid,
      day = day,
      age = age,
      treatment = ifelse(is.na(current_treatment), "none", current_treatment),
      state = state,
      duration_in_state = duration_in_state,
      ppfev = ppfev
    )
  }
  
  rbindlist(daily_data)
}

# Run Simulation ------------------------------------------------------------
cat("Simulating patient trajectories...\n")
simulation_data <- rbindlist(lapply(1:nrow(patients), function(i) {
  if (i %% 10 == 0) cat("  Patient", i, "of", n_patients, "\n")
  simulate_patient(patients[i])
}))

# Summary Statistics --------------------------------------------------------
cat("\n=== SIMULATION SUMMARY ===\n")
cat("Total patients:", n_patients, "\n")
cat("Total days:", n_days, "\n")

# Count exacerbations per patient
simulation_data[, state_change := c(FALSE, state[-1] != state[-.N]), by = patient_id]
exacerbations <- simulation_data[state_change == TRUE & state == "exacerbation", 
                                 .(n_exacerbations = .N), by = patient_id]

# Fill in patients with zero exacerbations
all_patients <- data.table(patient_id = 1:n_patients)
exacerbations <- merge(all_patients, exacerbations, all.x = TRUE)
exacerbations[is.na(n_exacerbations), n_exacerbations := 0]

cat("\nExacerbations per patient:\n")
print(summary(exacerbations$n_exacerbations))
cat("Mean exacerbations per patient:", mean(exacerbations$n_exacerbations), "\n")

# Calculate exacerbation episode durations
simulation_data[, episode := rleid(state), by = patient_id]
exacerb_episodes <- simulation_data[state == "exacerbation", 
                                    .(duration = .N, 
                                      mean_ppfev = mean(ppfev),
                                      treatment = first(treatment)), 
                                    by = .(patient_id, episode)]

if (nrow(exacerb_episodes) > 0) {
  cat("\nExacerbation episode durations (days):\n")
  print(summary(exacerb_episodes$duration))
  
  cat("\nTreatment distribution across exacerbations:\n")
  print(table(exacerb_episodes$treatment))
  cat("\nTotal exacerbation episodes:", nrow(exacerb_episodes), "\n")
  
  cat("\nMean exacerbation duration by treatment:\n")
  print(exacerb_episodes[, .(
    mean_duration = mean(duration),
    sd_duration = sd(duration),
    count = .N
  ), by = treatment])
}

# ppFEV summary
cat("\nppFEV summary:\n")
cat("Overall:\n")
print(summary(simulation_data$ppfev))

cat("\nBy state:\n")
print(simulation_data[, .(
  mean_ppfev = mean(ppfev),
  sd_ppfev = sd(ppfev),
  min_ppfev = min(ppfev),
  max_ppfev = max(ppfev)
), by = state])

# Save Data -----------------------------------------------------------------
fwrite(simulation_data, "patient_simulation_data.csv")
cat("\nData saved to: patient_simulation_data.csv\n")

# Create Visualizations -----------------------------------------------------
cat("\nCreating visualizations...\n")

# 1. Patient trajectories for first 5 patients
sample_patients <- unique(simulation_data$patient_id)[1:6]
plot_data <- simulation_data[patient_id %in% sample_patients]

p1 <- ggplot(plot_data, aes(x = day, y = ppfev)) +
  geom_point(aes(, color = state), size = 0.8) +
  geom_line(lwd = 0.3) +
  facet_wrap(~patient_id, ncol = 2,
             labeller = labeller(patient_id = function(x) paste("Patient", x))) +
  scale_color_manual(values = c("healthy" = "#2E7D32", "exacerbation" = "#C62828")) +
  labs(title = "Patient Trajectories: ppFEV Over Time",
       subtitle = "First 6 patients",
       x = "Day", y = "ppFEV (%)", color = "State") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"))

p1

# 2. Distribution of exacerbations
p2 <- ggplot(exacerbations, aes(x = n_exacerbations)) +
  geom_bar(fill = "#1976D2", alpha = 0.8) +
  labs(title = "Distribution of Exacerbations per Patient",
       subtitle = paste0("Mean = ", round(mean(exacerbations$n_exacerbations), 2)),
       x = "Number of Exacerbations", y = "Number of Patients") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14))

p2

# 3. Mean ppFEV by treatment over time (only during exacerbations)
exacerb_data <- simulation_data[state == "exacerbation" & treatment != "none"]
mean_ppfev <- exacerb_data[, .(mean_ppfev = mean(ppfev), n = .N), by = .(day, treatment)]

p3 <- ggplot(mean_ppfev, aes(x = day, y = mean_ppfev, color = treatment)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Mean ppFEV by Treatment Arm During Exacerbations",
       subtitle = "(Treatment randomized at each exacerbation)",
       x = "Day", y = "Mean ppFEV During Exacerbation (%)", color = "Treatment") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))

p3

# 4. Time in exacerbation
state_summary <- simulation_data[, .(
  days_healthy = sum(state == "healthy"),
  days_exacerbation = sum(state == "exacerbation")
), by = patient_id]
state_summary[, prop_exacerbation := days_exacerbation / (days_healthy + days_exacerbation)]

p4 <- ggplot(state_summary, aes(x = prop_exacerbation * 100)) +
  geom_histogram(fill = "#1976D2", alpha = 0.7, color = "black", bins = 20) +
  geom_vline(xintercept = mean(state_summary$prop_exacerbation) * 100, 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = mean(state_summary$prop_exacerbation) * 100 + 1, 
           y = Inf, vjust = 2,
           label = paste0("Mean = ", round(mean(state_summary$prop_exacerbation) * 100, 1), "%"),
           color = "red", fontface = "bold") +
  labs(title = "Distribution of Time Spent in Exacerbation State",
       x = "Proportion of Time in Exacerbation (%)", 
       y = "Number of Patients") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14))

p4

# 5. Exacerbation duration by treatment
if (nrow(exacerb_episodes) > 0) {
  exacerb_episodes_clean <- exacerb_episodes[treatment != "none"]
  
  p5 <- ggplot(exacerb_episodes_clean, aes(x = treatment, y = duration, fill = treatment)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("standard" = "#2E7D32", 
                                 "delayed" = "#F57C00", 
                                 "deferred" = "#C62828")) +
    labs(title = "Exacerbation Duration by Randomized Treatment",
         subtitle = "(Treatment assigned at exacerbation onset)",
         x = "Treatment", y = "Exacerbation Duration (days)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14)) +
    scale_x_discrete(labels = c("deferred" = "Deferred", 
                                "delayed" = "Delayed",
                                "standard" = "Standard"))
  
  p5
}

cat("\nAll visualizations saved!\n")
cat("\nSimulation complete!\n")

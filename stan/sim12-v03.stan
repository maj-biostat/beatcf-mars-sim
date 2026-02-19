// Weibull-PH with Gamma frailty — marginalised over frailty
data {
  int<lower=0> N_obs;
  int<lower=0> N_id;
  vector[N_obs] y_obs;
  array[N_obs] int<lower=1, upper=N_id> id;   // must be 1-indexed and contiguous
  int<lower=1> P;
  matrix[N_obs, P] X;
  real<lower=0> pri_s_u;
}
parameters {
  real<lower=0> shape;       // Weibull shape (alpha)
  real          b_0;         // intercept
  vector[P]     b;           // covariate effects
  real<lower=0> u_a;         // Gamma frailty shape (= rate, so mean = 1)
}
transformed parameters {
  // linear predictor per obs (no frailty here — that's marginalised)
  vector[N_obs] log_mu;
  for (i in 1:N_obs){
    log_mu[i] = b_0 + dot_product(X[i], b);  
  }
    
}
model {
  // --- priors ---
  target += exponential_lpdf(shape | 1);
  target += normal_lpdf(b_0   | 0, 3);
  target += normal_lpdf(b     | 0, 3);
  target += exponential_lpdf(u_a | pri_s_u);

  // --- marginalised gamma-frailty Weibull-PH likelihood ---
  // Accumulate per-subject sufficient statistics
  vector[N_id] sum_log_density;  // sum of log(alpha * mu_ij * y_ij^(alpha-1))
  vector[N_id] S;                // sum of mu_ij * y_ij^alpha  (cumulative hazard per subject)
  array[N_id]  int n_i;          // event count per subject

  sum_log_density = rep_vector(0.0, N_id);
  S               = rep_vector(0.0, N_id);
  for (k in 1:N_id) {
    n_i[k] = 0;
  }

  for (i in 1:N_obs) {
    int sid = id[i];
    real y_alpha = pow(y_obs[i], shape);
    real mu_i    = exp(log_mu[i]);

    sum_log_density[sid] += log(shape) + log_mu[i] + (shape - 1) * log(y_obs[i]);
    S[sid]               += mu_i * y_alpha;
    n_i[sid]             += 1;
  }

  // Per-subject marginal contribution
  for (k in 1:N_id) {
    real nk = n_i[k];
    target += sum_log_density[k]
              + lgamma(nk + u_a) - lgamma(u_a)
              + u_a * log(u_a)
              - (nk + u_a) * log(S[k] + u_a);
  }
}
generated quantities {
  // Posterior mean frailty per subject (E[u | data] = (n_i + u_a)/(S_i + u_a))
  // Can be computed here if needed for diagnostics
}

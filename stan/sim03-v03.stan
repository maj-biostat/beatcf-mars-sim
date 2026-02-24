data {
  // number of unique pt
  int<lower=1> N;        
  // number of timepoints
  int<lower=1> J;               
  // observed outcome (ppFEV1)
  matrix[N, J] y;               
  // age at baseline (years)
  vector[N] age;             
  // time in fraction of fu (0, 0.25, 0.5, 0.75, 1)
  vector[J] t_obs;      
  // trt allocation
  array[N] int trt;
}
transformed data{
  vector[N] ln_age;
  for(i in 1:N){
    ln_age[i] = log(age[i]);  
  }
}
parameters {
  real b0;
  real b_ln_age;
  real b_t;
  vector[2] b_trt_raw;                  
  
  real<lower=0> sigma;  
  // unconstrained AR(1) correlation      
  real rho_un; 
}
transformed parameters {
  vector[3] b_trt;
  b_trt[1] = 0.0;
  b_trt[2:3] = b_trt_raw;
  
  real<lower=-1,upper=1> rho = tanh(rho_un);
  
}
model {
  // Priors
  target += normal_lpdf(b0 | 100, 5);
  target += normal_lpdf(b_ln_age | 0, 5);
  target += normal_lpdf(b_t | 0, 5);
  target += normal_lpdf(b_trt | 0, 5);
  
  target += exponential_lpdf(sigma | 0.2);
  target += normal_lpdf(rho_un | 0, 1);

  // subject-specific intercept part
  vector[N] base;           
  for (i in 1:N) {
    base[i] = b0 + b_trt[trt[i]] + b_ln_age * ln_age[i];
  }
  // mean matrix
  // vectorized: adds scalar to each row
  matrix[N, J] mu;          
  for (j in 1:J) {
    mu[, j] = base + b_t * t_obs[j];  
  }
  // elementwise residuals
  matrix[N, J] resid = y - mu;   

  // first column: resid[,1] ~ Normal(0, sigma) (vectorized over N)
  {
    vector[N] r1 = to_vector(resid[, 1]);
    target += normal_lpdf(r1 | 0, sigma);
  }
    
  // innovations SD for conditional steps (given marginal parametrization)
  real innov_sd = sigma * sqrt(1 - square(rho));

  // for j = 2..J compute eta_j = resid[,j] - rho * resid[,j-1], vectorized over N
  for (j in 2:J) {
    vector[N] eta_j = to_vector(resid[, j]) - rho * to_vector(resid[, j - 1]);
    target += normal_lpdf(eta_j | 0, innov_sd);
  }
  
}



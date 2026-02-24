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
  
  matrix[J, J] R;
  for (i in 1:J){
    for (j in 1:J){
      R[i,j] = pow(rho, abs(i-j));  
    }
  }
    
  matrix[J, J] Sigma = sigma^2 * R;
  
}
model {
  // Priors
  target += normal_lpdf(b0 | 100, 10);
  target += normal_lpdf(b_ln_age | 0, 10);
  target += normal_lpdf(b_t | 0, 10);
  target += normal_lpdf(b_trt | 0, 10);
  
  target += exponential_lpdf(sigma | 0.2);
  target += normal_lpdf(rho_un | 0, 1);

  // Likelihood: multivariate normal
  for (i in 1:N) {
    vector[J] mu_i;
    // vectorized mean for subject i
    mu_i = rep_vector(b0 + b_trt[trt[i]] + b_ln_age * ln_age[i], J)
           + b_t * t_obs;
    
    target += multi_normal_lpdf(to_vector(y[i]') | mu_i, Sigma);
  }
    
  
}



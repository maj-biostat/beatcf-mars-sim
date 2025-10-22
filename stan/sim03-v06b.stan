data {
  // number of unique subjects
  int<lower=1> N;        
  // number of timepoints
  int<lower=1> J;               
  // observed outcome (ppFEV1)
  matrix[N,J] y;
  // age at baseline (years)
  vector[N] age;             
  // time in fraction of follow-up (0, 0.25, 0.5, 0.75, 1)
  vector[J] t_obs;      
  // treatment allocation
  array[N] int trt;
}
transformed data {
  vector[N] ln_age;
  for (i in 1:N) {
    ln_age[i] = log(age[i]);
  }
}
parameters {
  real b0;
  real b_age;
  real b_t;
  vector[2] b_trt_raw;

  // marginal SD - see distinction in notes
  real<lower=0> sigma;    
  // unconstrained AR(1) correlation parameter
  real rho_un;            

}
transformed parameters {
  vector[3] b_trt;
  b_trt[1] = 0.0;
  b_trt[2:3] = b_trt_raw;

  real<lower=-1,upper=1> rho = tanh(rho_un);

  // // Build AR(1) correlation matrix
  // matrix[J, J] R;
  // for (r in 1:J) {
  //   for (c in 1:J) {
  //     R[r, c] = pow(rho, abs(r - c));
  //   }
  // }

}

model {
  // Priors
  target += normal_lpdf(b0 | 100, 5);
  target += normal_lpdf(b_age | 0, 5);
  target += normal_lpdf(b_t | 0, 5);
  target += normal_lpdf(b_trt | 0, 5);

  target += exponential_lpdf(sigma | 0.5);
  target += normal_lpdf(rho_un | 0, 1);
  
  matrix[N,J] mu;        
  real innov_sd = sigma * sqrt(1 - (rho*rho));
  
  for(j in 1:J){
    mu[,j] = b0 +
      // constant trt effect
      b_trt[trt] +
      // decline based on log baseline age
      b_age * ln_age  +
      // linear time trend
      b_t * t_obs[j];  
  }
  
  matrix[N,J] resid = y - mu;
  
  target += normal_lpdf(y[, 1] | mu[, 1], sigma);
  
  for(j in 2:J){
    target += normal_lpdf(y[, j] | mu[, j] + rho * resid[, j-1], innov_sd);
  }
  
}

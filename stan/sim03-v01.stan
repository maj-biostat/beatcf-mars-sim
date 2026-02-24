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
  // degrees of freedom for multivariate t (nu > 2)
  real<lower=2> nu;                     
  
  real<lower=0> sigma;             
  real<lower=-1,upper=1> rho;  // AR(1) correlation
  
  // subject-level scale parameters (gamma latent)
  vector<lower=0>[N] lambda;            
}
transformed parameters {
  
  matrix[J, J] Sigma;
  for (i in 1:J){
    for (j in 1:J){
      Sigma[i,j] = pow(rho, abs(i-j));  
    }
  }
    
  matrix[J, J] L = cholesky_decompose(Sigma) * sigma;
  
  vector[3] b_trt;
  b_trt[1] = 0.0;
  b_trt[2:3] = b_trt_raw;
}
model {
  // Priors
  target += normal_lpdf(b0 | 100, 10);
  target += normal_lpdf(b_ln_age | 0, 10);
  target += normal_lpdf(b_t | 0, 10);
  target += normal_lpdf(b_trt | 0, 10);
  
  target += exponential_lpdf(sigma | 0.2);
  target += exponential_lpdf(nu - 2 | 1.0 / 30.0); 
  
  target += uniform_lpdf(rho | -1, 1);

  // Prior for lambda implied by mixture representation:
  // lambda_i ~ Gamma(nu/2, nu/2)
  for (i in 1:N){
    lambda[i] ~ gamma(nu / 2.0, nu / 2.0);  
  }
    
  // Likelihood via scale mixture:
  for (i in 1:N) {
    vector[J] mu_i;
    vector[J] yi = to_vector(y[i]'); 
    
    for (j in 1:J){
      mu_i[j] = b0 + b_ln_age * ln_age[i] + b_t * t_obs[j] + b_trt[trt[i]];  
    }
    
    // scaled cholesky for conditional multivariate normal
    yi ~ multi_normal_cholesky(mu_i, L / sqrt(lambda[i]));
  }
}



data {
  // number of observations (repeat measures on pts at 3, 6, 9, 12 months)
  // all pt must have at least one follow up 
  // (baseline only pt dropped from analysis)
  int<lower=1> N;               
  // observed outcome (ppFEV1)
  vector[N] y;
  // design excl intercept
  
  // y0 baseline ppFEV1 (scaled)
  // log(age0)
  // t_obs0.5
  // t_obs0.75
  // t_obs1
  // trt2
  // trt3
  // t_obs0.5:trt2
  // t_obs0.75:trt2
  // t_obs1:trt2
  // t_obs0.5:trt3
  // t_obs0.75:trt3
  // t_obs1:trt3
  
  int<lower=0> P;
  matrix[N, P] X;
  
  // ix for fu1 
  int N_1;
  array[N_1] int ix_1;
  int N_2;
  array[N_2] int ix_2;
  int N_3;
  array[N_3] int ix_3;
  int N_4;
  array[N_4] int ix_4;
}
transformed data {
  array[N_2] int ix_2_prev;
  array[N_3] int ix_3_prev;
  array[N_4] int ix_4_prev;
  
  for(i in 1:N_2) ix_2_prev[i] = ix_2[i] - 1;
  for(i in 1:N_3) ix_3_prev[i] = ix_3[i] - 1;
  for(i in 1:N_4) ix_4_prev[i] = ix_4[i] - 1;
}
parameters {
  real b_0;
  vector[P] b;

  // marginal SD - see distinction in notes
  real<lower=0> sigma;    
  // unconstrained AR(1) correlation parameter
  real rho_un;            

}
transformed parameters {
  real<lower=-1,upper=1> rho = tanh(rho_un);
}
model {
  // Priors
  target += normal_lpdf(b_0 | 100, 5);
  target += normal_lpdf(b | 0, 5);

  target += exponential_lpdf(sigma | 0.5);
  target += normal_lpdf(rho_un | 0, 1);
  
  vector[N] mu = b_0 + X * b;        
  real innov_sd = sigma * sqrt(1 - (rho*rho));
  
  vector[N] resid = y - mu;
  
  target += normal_lpdf(y[ix_1] | mu[ix_1], sigma);

  target += normal_lpdf(y[ix_2] | mu[ix_2] + rho * resid[ix_2_prev], innov_sd);
  target += normal_lpdf(y[ix_3] | mu[ix_3] + rho * resid[ix_3_prev], innov_sd);
  target += normal_lpdf(y[ix_4] | mu[ix_4] + rho * resid[ix_4_prev], innov_sd);
  
  
}
generated quantities{
  
  
}

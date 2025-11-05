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
  int N_5;
  array[N_5] int ix_5;
  int N_6;
  array[N_6] int ix_6;
  int N_7;
  array[N_7] int ix_7;
  int N_8;
  array[N_8] int ix_8;
  int N_9;
  array[N_9] int ix_9;
  int N_10;
  array[N_10] int ix_10;
  int N_11;
  array[N_11] int ix_11;
  int N_12;
  array[N_12] int ix_12;
}
transformed data {
  array[N_2] int ix_2_prev;
  array[N_3] int ix_3_prev;
  array[N_4] int ix_4_prev;
  array[N_5] int ix_5_prev;
  array[N_6] int ix_6_prev;
  array[N_7] int ix_7_prev;
  array[N_8] int ix_8_prev;
  array[N_9] int ix_9_prev;
  array[N_10] int ix_10_prev;
  array[N_11] int ix_11_prev;
  array[N_12] int ix_12_prev;
  
  for(i in 1:N_2) ix_2_prev[i] = ix_2[i] - 1;
  for(i in 1:N_3) ix_3_prev[i] = ix_3[i] - 1;
  for(i in 1:N_4) ix_4_prev[i] = ix_4[i] - 1;
  for(i in 1:N_5) ix_5_prev[i] = ix_5[i] - 1;
  for(i in 1:N_6) ix_6_prev[i] = ix_6[i] - 1;
  for(i in 1:N_7) ix_7_prev[i] = ix_7[i] - 1;
  for(i in 1:N_8) ix_8_prev[i] = ix_8[i] - 1;
  for(i in 1:N_9) ix_9_prev[i] = ix_9[i] - 1;
  for(i in 1:N_10) ix_10_prev[i] = ix_10[i] - 1;
  for(i in 1:N_11) ix_11_prev[i] = ix_11[i] - 1;
  for(i in 1:N_12) ix_12_prev[i] = ix_12[i] - 1;
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
  target += normal_lpdf(y[ix_5] | mu[ix_5] + rho * resid[ix_5_prev], innov_sd);
  target += normal_lpdf(y[ix_6] | mu[ix_6] + rho * resid[ix_6_prev], innov_sd);
  target += normal_lpdf(y[ix_7] | mu[ix_7] + rho * resid[ix_7_prev], innov_sd);
  target += normal_lpdf(y[ix_8] | mu[ix_8] + rho * resid[ix_8_prev], innov_sd);
  target += normal_lpdf(y[ix_9] | mu[ix_9] + rho * resid[ix_9_prev], innov_sd);
  target += normal_lpdf(y[ix_10] | mu[ix_10] + rho * resid[ix_10_prev], innov_sd);
  target += normal_lpdf(y[ix_11] | mu[ix_11] + rho * resid[ix_11_prev], innov_sd);
  target += normal_lpdf(y[ix_12] | mu[ix_12] + rho * resid[ix_12_prev], innov_sd);
  
  
}
generated quantities{
  
  
}

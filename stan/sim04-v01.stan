data {
  // number of observations
  int<lower=1> N;               
  
  // observed outcome (ppFEV1)
  vector[N] y;
  
  // number of subjects
  int<lower=1> N_S;
  array[N] int id;
  vector[N] t;
  
  // design excl intercept
  // time
  // log baseline age
  // trt2 indicator
  // trt3 indicator
  
  int<lower=1> P;  // number of population-level effects
  matrix[N, P] X;  // population-level design matrix
  
  
}
transformed data {
  
  vector[P] mu_X;
  for(i in 1:P){
    mu_X[i] = mean(X[, i]);
  }
  matrix[N, P] X_c;
  for(i in 1:P){
    X_c[, i] = X[, i] - mu_X[i];
  }
}
parameters {
  real b_0_tmp;
  vector[P] b;

  vector[N_S] z0;
  vector[N_S] z1;
  
  // random intercept 
  real<lower=0> s_u0;
  // random time trend
  real<lower=0> s_u1;
           
  // residual
  real<lower=0> s_e;

}
transformed parameters {
  vector[N_S] u0;
  vector[N_S] u1;
  
  u0 = z0 * s_u0;
  u1 = z1 * s_u1;
}
model {
  // Priors
  target += normal_lpdf(b_0_tmp | 80, 10);
  target += normal_lpdf(b | 0, 5);

  target += exponential_lpdf(s_u0 | 1);
  target += exponential_lpdf(s_u1 | 1);
  
  target += exponential_lpdf(s_e | 1);
  
  target += normal_lpdf(z0 | 0, 1);
  target += normal_lpdf(z1 | 0, 1);
  
  vector[N] mu = b_0_tmp + X_c * b + u0[id];
  for(i in 1:N){
    mu[i] = mu[i] + u1[id[i]] * t[i];
  }
  
  target += normal_lpdf(y | mu, s_e);
  
}
generated quantities{
  real b_0 = b_0_tmp - dot_product(mu_X, b);
}

data{
  int N;
  vector[N] y;
  
  int Kx;
  // y_pre and trt indicators
  matrix[N, Kx] X;
  
  // priors
  // mean, sd
  vector[2] pri_b_0;
  // means for each beta
  vector[Kx] pri_b_mu;
  // sd for each beta
  vector[Kx] pri_b_sd;
  // rho for exponential
  real pri_se;
  
  int prior_only;
}
transformed data {
  matrix[N, Kx] Xc;  
  vector[Kx] mu_X;
  for (i in 1:Kx) {
    mu_X[i] = mean(X[, i]);
    Xc[, i] = X[, i] - mu_X[i];
  }
}
parameters{
  real b_0_tmp;
  vector[Kx] b;
  real<lower=0> se;
}
transformed parameters{
  
  vector[N] mu;
  mu = b_0_tmp + Xc * b;
  
}
model{
  
  target += normal_lpdf(b_0_tmp | pri_b_0[1], pri_b_0[2]);
  for(i in 1:Kx){
    target += normal_lpdf(b | pri_b_mu[i], pri_b_sd[i]);  
  }
  target += exponential_lpdf(se | pri_se);
  
  if(!prior_only){
    target += normal_lpdf(y | mu, se);  
  }
  
}
generated quantities{
  
  real b_0 = b_0_tmp - dot_product(mu_X, b);

  vector[N] w = dirichlet_rng(rep_vector(1, N));

  // In this setting, b_trt[2] and b_trt[3] will be identical to
  // delta_2_1 and delta_3_1.
  vector[N] eta_1 = b_0 + b[1]*X[, 1] ;
  vector[N] eta_2 = b_0 + b[1]*X[, 1] + b[2];
  vector[N] eta_3 = b_0 + b[1]*X[, 1] + b[3];

  real mu_1;
  real mu_2;
  real mu_3;
  
  mu_1 = w' * eta_1   ;
  mu_2 = w' * eta_2   ;
  mu_3 = w' * eta_3   ;

  // start with the marginal effect of treatment
  real delta_2_1 = mu_2 - mu_1;
  real delta_3_1 = mu_3 - mu_1;
  real delta_3_2 = mu_3 - mu_2;
  
}

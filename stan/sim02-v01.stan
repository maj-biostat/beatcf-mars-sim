data{
  int N;
  vector[N] y;
  vector[N] y_pre;
  
  // intervention
  array[N] int trt;
  
  // priors
  vector[2] pri_b_0;
  vector[2] pri_b_trt;
  vector[2] pri_b_pre;
  
  int prior_only;
}
parameters{
  real b_0;
  vector[2] b_trt_raw;
  real<lower=0> se;
  real b_pre;
}
transformed parameters{
  
  vector[3] b_trt;
  vector[N] eta;
  
  b_trt[1] = 0.0;
  
  b_trt[2:3] = b_trt_raw;
  
  // in practice would assume a non linear effect on pre
  eta = b_0 + b_pre*y_pre + b_trt[trt];
}
model{
  
  target += normal_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  target += normal_lpdf(b_trt_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(b_trt_raw | pri_b_pre[1], pri_b_pre[2]);
  
  if(!prior_only){
    target += normal_lpdf(y | eta, se);  
  }
  
}
generated quantities{
  
  vector[N] w = dirichlet_rng(rep_vector(1, N));

  vector[N] eta_1 = b_0 + b_pre*y_pre + b_trt[1];
  vector[N] eta_2 = b_0 + b_pre*y_pre + b_trt[2];
  vector[N] eta_3 = b_0 + b_pre*y_pre + b_trt[3];

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

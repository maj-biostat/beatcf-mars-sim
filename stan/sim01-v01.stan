data{
  int N;
  array[N] int y;
  array[N] int n;
  
  // intervention
  array[N] int trt;
  
  // lung function (< 70% pFEV1) = 1, (>= 70% pFEV1) = 2
  array[N] int lung;
  // colonisation (no pseudo) = 1, (pseudo) = 2 
  array[N] int pseudo;
  
  // priors
  vector[2] pri_a;
  vector[2] pri_b_trt;
  vector[2] pri_b_lung;
  vector[2] pri_b_pseudo;
  
  int prior_only;
}
parameters{
  real a;
  vector[2] b_trt_raw;
  real b_lung_raw;
  real b_pseudo_raw;
  
}
transformed parameters{
  
  vector[3] b_trt;
  vector[2] b_lung;
  vector[2] b_pseudo;
  vector[N] eta;
  
  b_trt[1] = 0.0;
  b_lung[1] = 0.0;
  b_pseudo[1] = 0.0;
  
  b_trt[2:3] = b_trt_raw;
  b_lung[2] = b_lung_raw;
  b_pseudo[2] = b_pseudo_raw;
  
  eta = a + b_trt[trt] + b_lung[lung] + b_pseudo[pseudo];
  
}
model{
  
  target += logistic_lpdf(a | pri_a[1], pri_a[2]);
  target += normal_lpdf(b_trt_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(b_lung_raw | pri_b_lung[1], pri_b_lung[2]);
  target += normal_lpdf(b_pseudo_raw | pri_b_pseudo[1], pri_b_pseudo[2]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }
  
}
generated quantities{
  
  vector[N] w = dirichlet_rng(to_vector(n));

  vector[N] mu_1 = inv_logit(a + b_trt[1] + b_lung[lung] + b_pseudo[pseudo]);
  vector[N] mu_2 = inv_logit(a + b_trt[2] + b_lung[lung] + b_pseudo[pseudo]);
  vector[N] mu_3 = inv_logit(a + b_trt[3] + b_lung[lung] + b_pseudo[pseudo]);

  real p_1 = w' * mu_1   ;
  real p_2 = w' * mu_2   ;
  real p_3 = w' * mu_3   ;

  // start with the marginal effect of treatment
  real rd_2_1 = p_2 - p_1;
  real rd_3_1 = p_3 - p_1;
  real rd_3_2 = p_3 - p_2;
  
}

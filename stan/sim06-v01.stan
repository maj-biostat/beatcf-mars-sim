data {
  int N;
  int N_pt;
  
  // number of periods
  int P; 
  array[N] int period;
  array[N] int evt;
  
  int K;
  array[N] int trt;
  
  // data needs to be rejigged so that ids are 1:N_pt
  array[N] int id;
}

parameters {
  vector[P] alpha;
  real<lower = 0> s_u;
  vector[N_pt] u_raw;
  vector[K-1] b_trt_raw;
}

transformed parameters {
  vector[N_pt] u = u_raw * s_u;
  vector[K] b_trt;
  
  b_trt[1] = 0.0;
  b_trt[2:K] = b_trt_raw;
  
}

model {
  
  target += normal_lpdf(u_raw | 0, 1);
  target += exponential_lpdf(s_u | 1);
  target += normal_lpdf(b_trt_raw | 0, 3);
  
  target += bernoulli_lpmf(evt | inv_logit(alpha[period] + b_trt[trt] + u[id]));
  
}

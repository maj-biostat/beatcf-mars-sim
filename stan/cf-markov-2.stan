data {
  int<lower=1> N;                 // observations
  int<lower=1> P; // matrix cols
  
  matrix[N, P] X;
  array[N] int y; // 1 none, 2 mild, 3 severe

  int ix_trt_2;
  int  ix_trt_3 ;
  int  ix_prev_2;
  int  ix_prev_3;
  int  ix_time_1;
  int  ix_time_2;
  int  ix_gap ;
  int  ix_trt_time_2;
  int  ix_trt_time_3;
  
  real mu_days;
  real sd_days;
  
  int<lower=1> T_pred; // usually 28
}
transformed data{
  vector[1] zero1 = rep_vector(0, 1);
}
parameters {
  ordered[2] alpha;
  vector[P] b;
}
transformed parameters {
}
model {

  // priors
  alpha ~ normal(0, 2);
  b ~ normal(0, 1);

  // likelihood
  
  vector[N] eta = X * b;
  y ~ ordered_logistic(eta, alpha);

}
generated quantities{
  
  // vector[P] g;
  // vector[2] alpha;
  // 
  // g[ix_no_trans] = b[ix_no_trans];
  // g[ix_time_1] = (b[ix_time_1] / sd_days) - (2.0 * mu_days * b[ix_time_2]) / pow(sd_days, 2);
  // g[ix_time_2] = b[ix_time_2] / pow(sd_days, 2);
  // 
  // // 4. Back-transform Interaction Coefficients (Rescaled by sigma)
  // g[ix_trt_time_1] = b[ix_trt_time_1] / sd_days;
  // g[ix_trt_time_2] = b[ix_trt_time_2] / sd_days;
  // 
  // // 5. Back-transform Treatment Main Effects (Shifted by the interaction)
  // g[ix_trt_1] = b[ix_trt_1] - (b[ix_trt_time_1] * mu_days) / sd_days;
  // g[ix_trt_2] = b[ix_trt_2] - (b[ix_trt_time_2] * mu_days) / sd_days;
  // 
  // 
  // for (k in 1:2) {
  //   alpha[k] = a[k] - ((b[ix_time_2] * pow(mu_days,2)) / pow(sd_days, 2) - (b[ix_time_1] * mu_days) / sd_days);
  // }
  
  // see sim18-back-transform-note.qmd
  vector[2] a;
  for (k in 1:2) {
    a[k] = alpha[k] - ((b[ix_time_2] * pow(mu_days,2)) / pow(sd_days, 2) - (b[ix_time_1] * mu_days) / sd_days);
  }

  vector[3] b_trt;
  b_trt[1] = 0.0;
  b_trt[2] = b[ix_trt_2] - (b[ix_trt_time_2] * mu_days) / sd_days;
  b_trt[3] = b[ix_trt_3] - (b[ix_trt_time_3] * mu_days) / sd_days;
  
  vector[3] b_prev;
  b_prev[1] = 0.0;
  b_prev[2] = b[ix_prev_2];
  b_prev[3] = b[ix_prev_3];
  
  real b_time_1 = (b[ix_time_1] / sd_days) - (2.0 * mu_days * b[ix_time_2]) / pow(sd_days, 2);
  real b_time_2 = b[ix_time_2] / pow(sd_days, 2);
  
  real b_gap = b[ix_gap];
  
  vector[3] b_trt_time;
  b_trt_time[1] = 0.0;
  b_trt_time[2] = b[ix_trt_time_2] / sd_days;
  b_trt_time[3] = b[ix_trt_time_3] / sd_days;
  
}

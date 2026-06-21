data {
  int<lower=1> N;                 // observations
  
  vector[N] time;
  array[N] int gap_time;  // 1 = 1 day, 2 = 7 day 
  
  array[N] int<lower=1, upper=3> trt; 
  
  array[N] int prev; // 1 none, 2 mild, 3 severe
  array[N] int y; // 1 none, 2 mild, 3 severe

  int<lower=1> T_pred; // usually 28
}
transformed data{
  vector[1] zero1 = rep_vector(0, 1);
}
parameters {

  ordered[2] a;

  vector[2] b_trt_raw;
  vector[2] b_prev_raw;
  real b_time_1;
  real b_time_2;
  vector[2] b_trt_time_raw;
  vector[1] b_gap_raw;

}

transformed parameters {

  vector[3] b_trt = append_row(zero1, b_trt_raw); 
  
  // treatment by time interaction as quadratic
  vector[3] b_trt_time = append_row(zero1, b_trt_time_raw); 
  
  vector[3] b_prev = append_row(zero1, b_prev_raw); 
  vector[2] b_gap = append_row(zero1, b_gap_raw); 
  
}

model {

  // priors
  a ~ normal(0, 2);

  b_trt_raw ~ normal(0, 1);
  b_prev_raw ~ normal(0, 1);
  b_time_1 ~ normal(0, 1);
  b_time_2 ~ normal(0, 1);
  b_trt_time_raw ~ normal(0, 1);
  b_gap_raw ~ normal(0, 1);

  // likelihood
  
  vector[N] eta = b_trt[trt] + b_prev[prev] + 
    b_time_1*time + b_time_2*pow(time, 2) +
    b_gap[gap_time] +
    b_trt_time[trt] .* time;
    
  y ~ ordered_logistic(eta, a);

}


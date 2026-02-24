// weibull ph model
functions {
  real weibull2_lpdf(vector y, real alpha, vector lambda) {
    int nn = num_elements(y);
    
    vector[nn] prob;
    vector[nn] lprob;
    
    for(i in 1:nn){
      // pdf
      prob[i] = alpha * lambda[i] * y[i]^(alpha - 1) * exp(- lambda[i]*y[i]^alpha);
      if(prob[i] < machine_precision()) {
        prob[i] = machine_precision();
      }
      lprob[i] = log(prob[i]);
    }
    
    return sum(lprob);
  }
  
  real weibull2_lcdf(vector y, real alpha, vector lambda) {
    int nn = num_elements(y);
    
    vector[nn] prob;
    vector[nn] lprob;
    
    // cdf
    for(i in 1:nn){
      // pdf
      prob[i] = 1-exp(-lambda[i] * y[i]^alpha);
      if(prob[i] < machine_precision()) {
        prob[i] = machine_precision();
      }
      lprob[i] = log(prob[i]);
    }
    
    return sum(lprob);
  }

  real weibull2_lccdf(vector y, real alpha, vector lambda) {
    int nn = num_elements(y);
    
    vector[nn] prob;
    vector[nn] lprob;
    // 1 - cdf = S
    
    for(i in 1:nn){
      // pdf
      prob[i] = exp(-lambda[i] * y[i]^alpha);
      if(prob[i] < machine_precision()) {
        prob[i] = machine_precision();
      }
      lprob[i] = log(prob[i]);
    }
    
    return sum(lprob);
  }
  
  real weibull2_rng(real alpha, real lambda) {
    real u;
    real z;
    u = uniform_rng(0, 1);
    z = (-log(u)/lambda)^(1/alpha);
    return z;
  }

  real weibull2_lb_rng(real alpha, real lambda, real lb, int dbg) {
    
    real p = 1-exp(-lambda * lb^alpha);
    real u;
    real z;

    if(is_nan(p)) p = 0.9999;
    else if (p >= 1) p = 0.9999;
    else if (p <= 0) p = 0.0;

    u = uniform_rng(p, 1);
    z = (-log(1-u)/lambda)^(1/alpha);

    return z;
  }
}
data {
  int<lower=0> N_obs;
  int N_id;
  vector[N_obs] y_obs;
  array[N_obs] int id;
  // ppfev (mean centred)
  int P;
  matrix[N_obs, P] X;
  
  real<lower=0> pri_s_u;
}
parameters {

  real<lower=0> shape;
  real b_0;
  
  vector[P] b;
  vector<lower=0>[N_id] u_id;
  
  real<lower=0> u_a;
}
transformed parameters {
  vector[N_obs] scale;
  
  // the mean of the frailty (u_a / u_r) needs to be constrained to 1 otherwise
  // shifts will get absorbed into b_0 and vice versa
  real<lower=0> u_r = u_a;
  
  // with exp(b_0) * exp(b_trt) = exp(a + b)
  for(i in 1:N_obs){
    scale[i] = u_id[id[i]] * exp(b_0 + X[i]*b);   
  }
  
  
}
model {
  target += exponential_lpdf(shape | 1);
  target += normal_lpdf(b_0 | 0, 3);
  target += normal_lpdf(b | 0, 3);
  
  target += exponential_lpdf(u_a | 1);
  target += gamma_lpdf(u_id | u_a, u_r);
  
  // observed event time
  target += weibull2_lpdf(y_obs | shape, scale);
  
}
generated quantities {
}


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
    
    real p = exp(1-exp(-lambda * lb^alpha));
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
  int<lower=0> N_obs_he;
  int N_id_he;
  vector[N_obs_he] y_obs_he;
  array[N_obs_he] int id_he;
  // ppfev (mean centred)
  vector[N_obs_he] ppfev_he;
  
  int<lower=0> N_obs_eh;
  int N_id_eh;
  vector[N_obs_eh] y_obs_eh;
  array[N_obs_eh] int id_eh;
  // ppfev (mean centred)
  vector[N_obs_eh] ppfev_eh;
  array[N_obs_eh] int trt_ix_eh;
  
}
parameters {
  // healthy -> exacerbation
  real<lower=0> shape_he;
  real b_he_0;
  real b_he_ppfev;
  vector[N_id_he] z_he;
  real<lower=0> s_he;
  
  // exacerbation -> healthy
  real<lower=0> shape_eh;
  real b_eh_0;
  real b_eh_ppfev;
  vector[2] b_eh_trt_z;
  vector[N_id_eh] z_eh;
  real<lower=0> s_eh;
}
transformed parameters {
  vector[N_obs_he] scale_he;
  vector[N_obs_he] u_he;
  
  vector[N_obs_eh] scale_eh;
  vector[3] b_eh_trt;
  vector[N_obs_eh] u_eh;
  
  u_he = s_he * z_he[id_he];
  
  b_eh_trt[1] = 0.0;
  b_eh_trt[2:3] = b_eh_trt_z;
  u_eh = s_eh * z_eh[id_eh];
  
  // with exp(b_0) * exp(b_trt) = exp(a + b)
  scale_he = exp(b_he_0 + b_he_ppfev * ppfev_he + u_he);
  
  scale_eh = exp(b_eh_0 + b_eh_ppfev * ppfev_eh + b_eh_trt[trt_ix_eh] + u_eh);
  
}
model {
  target += exponential_lpdf(shape_he | 1);
  target += normal_lpdf(b_he_0 | 0, 3);
  target += normal_lpdf(b_he_ppfev | 0, 3);
  target += exponential_lpdf(s_he | 3);
  target += normal_lpdf(z_he | 0, 1);
  
  target += exponential_lpdf(shape_eh | 1);
  target += normal_lpdf(b_eh_0 | 0, 3);
  target += normal_lpdf(b_eh_ppfev | 0, 3);
  target += normal_lpdf(b_eh_trt_z | 0, 3);
  target += exponential_lpdf(s_eh | 3);
  target += normal_lpdf(z_eh | 0, 1);
  
  // observed event time
  target += weibull2_lpdf(y_obs_he | shape_he, scale_he);
  target += weibull2_lpdf(y_obs_eh | shape_eh, scale_eh);
  
}
generated quantities {
  // array[N] real yrep;
  // 
  // for(n in 1:N) {
  //   if(v[n] == 0){
  //     yrep[n] = weibull2_lb_rng(alpha, lp[n], y[n], dbg);
  //   } else if(v[n] == 1){
  //     yrep[n] = weibull2_rng(alpha, lp[n]);
  //   }
  //   
  // }
}


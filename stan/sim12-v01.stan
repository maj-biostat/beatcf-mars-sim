// weibull ph model
functions {
  real weibull2_lpdf(real y, real alpha, real lambda) {
    real prob;
    real lprob;
    // pdf
    prob = alpha * lambda * y^(alpha - 1) * exp(- lambda*y^alpha);
    // if(dbg) print("weibull2_lpdf : ", prob);
    if(prob < machine_precision()) {
      prob = machine_precision();
    }
    lprob = log(prob);
    return lprob;
  }
  
  real weibull2_lcdf(real y, real alpha, real lambda) {
    real prob;
    real lprob;
    // cdf
    prob = 1-exp(-lambda * y^alpha);
    // if(dbg) print("weibull2_lcdf : ", prob);
    if(prob < machine_precision()) {
      prob = machine_precision();
    }
    lprob = log(prob);
    return lprob;
  }

  real weibull2_lccdf(real y, real alpha, real lambda) {
    real prob;
    real lprob;
    // 1 - cdf = S
    prob = exp(-lambda * y^alpha);
    // if(dbg) print("weibull2_lccdf : ", prob, " ", alpha, " ", lambda);
    if(prob < machine_precision()) {
      prob = machine_precision();
    }
    lprob = log(prob);
    return lprob;
  }
  
  // real weibull2_rng(real alpha, real lambda) {
  //   real u;
  //   real z;
  //   u = uniform_rng(0, 1);
  //   z = (-log(u)/lambda)^(1/alpha);
  //   return z;
  // }
  // 
  // real weibull2_lb_rng(real alpha, real lambda, real lb, int dbg) {
  //   real p = exp(weibull2_lcdf(lb | alpha, lambda, dbg));
  //   real u;
  //   real z;
  //   
  //   if(is_nan(p)) p = 0.9999;
  //   else if (p >= 1) p = 0.9999;
  //   else if (p <= 0) p = 0.0;
  // 
  //   u = uniform_rng(p, 1);      
  //   z = (-log(1-u)/lambda)^(1/alpha);
  //   
  //   return z;
  // }
}
data {
  int<lower=0> N_obs;
  // int<lower=0> N_cens;
  array[N_obs] real y_obs;
  // array[N_cens] real y_cens;
  
  vector[N_obs] trt_obs;
  // array[N_cens] trt_cens;
  
}
parameters {
  real<lower=0> shape;
  real b_0;
  real b_trt;
}
transformed parameters {
  // real lambda0;
  vector[N_obs] scale_obs;
  
  # exp(b_0) * exp(b_trt)
  scale_obs = exp(b_0 + b_trt * trt_obs);
  // scale_cens = exp(X*beta);
  
}
model {
  target += exponential_lpdf(shape | 1);
  target += normal_lpdf(b_0 | 0, 3);
  target += normal_lpdf(b_trt | 0, 3);
  
  // observed event time
  for(i in 1:N_obs){
    target += weibull2_lpdf(y_obs[i] | shape, scale_obs[i]);  
  }
  
  
  // right-censored
  // target += weibull2_lccdf(y_cens | shape, scale_cens);
  
  
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


data {
  // number of observations
  int<lower=1> N;
  // total number of unique participants across both healthy and exacerbation
  // states
  int<lower=1> N_id;
  array[N] int id;
  
  array[N] int<lower=0,upper=1> y;
  
  // 1 represents health, 2 represents exacerbation
  array[N] int<lower=1,upper=2> state;
  
  // bin associated with current day of state   
  array[N] int<lower=1> bin;

  // number of bins in 
  int<lower=1> N_he_bin;
  int<lower=1> N_eh_bin;
  
  // indicator for first healthy period
  array[N] int i_entry;
  // design matrix
  int<lower=1> P;
  // ppfev at baseline, then treatment options
  matrix[N, P] X;
  int trt_defer_col;
  int trt_discont_col;
}

parameters {
  // piecewise baseline
  vector[N_he_bin] a_he;
  vector[N_eh_bin] a_eh;
  
  vector[P] b_he;
  vector[P] b_eh;
  
  vector[N_id] z_he;
  vector[N_id] z_eh;
  
  real<lower=0> sd_he;
  real<lower=0> sd_eh;
}

model {
  
  target += normal_lpdf(a_he | 0, 2);
  target += normal_lpdf(a_eh | 0, 2);
  
  target += normal_lpdf(b_he | 0, 1);
  target += normal_lpdf(b_eh | 0, 1);
  
  target += normal_lpdf(z_he | 0, 1);
  target += normal_lpdf(z_eh | 0, 1);
  
  target += exponential_lpdf(sd_he | 1);
  target += exponential_lpdf(sd_eh | 1);
  
  for (i in 1:N) {
    real eta;
    
    if(state[i]==1){
      
      if(i_entry[i] == 1){
        // exclude carry over effects
        eta = a_he[bin[i]] + X[i, 1]  * b_he[1] + z_he[id[i]] * sd_he;
      } else {
        eta = a_he[bin[i]] + X[i, ]  * b_he + z_he[id[i]] * sd_he;
      }
      
    } else{
      eta = a_eh[bin[i]] + X[i,]  * b_eh + z_eh[id[i]] * sd_eh;
    }
    
    target += bernoulli_logit_lpmf(y[i] | eta);
    
  }
}



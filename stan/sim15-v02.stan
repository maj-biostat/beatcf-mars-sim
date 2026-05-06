// revised to run based on start stop data
data {
  // number of observations (segments)
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
  vector[N] len_seg;
  // design matrix
  int<lower=1> P;
  // ppfev at baseline, then treatment options
  matrix[N, P] X;
  int trt_defer_col;
  int trt_discont_col;
  
  real pri_sd_he;
  real pri_sd_eh;
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
  
  target += exponential_lpdf(sd_he | 1/pri_sd_he);
  target += exponential_lpdf(sd_eh | 1/pri_sd_he);
  
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
    
    real p = inv_logit(eta);
    
    if(y[i] == 1){
      // here we have several daily interval contributions where no event 
      // occurred, followed by an event in the final interval, e.g.
      // (0, 1], (1, 2], (2, 3] no event
      // (3, 4] event (assumed to occur at end of interval)
      target += (len_seg[i] - 1) * log1m(p) + log(p);
    } else {
      // here we just have several daily intervals of no events.
      target += len_seg[i] * log1m(p);
    }
    
    
  }
}



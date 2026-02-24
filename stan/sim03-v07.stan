data {
  // number of rows (for complete data is Q x J)
  int<lower=1> N;  
  
  // number of unique pt
  int<lower=1> Q;
  // start and end indices for each pts obs
  // first obs is not modelled, we just condition on it
  array[Q] int s_id;
  array[Q] int e_id;
  
  // observed outcome (ppFEV1)
  vector[N] y;               
  
  // number of cov excl intercept
  int<lower = 0> P;
  matrix[N, P] X;
  
}
parameters {
  real b0;
  
  vector[P] b;
  
  real<lower=0> sigma;    
  // unconstrained autocorrel
  real rho_un;          

}
transformed parameters {
  
  real<lower=-1,upper=1> rho = tanh(rho_un);

}

model {
  // Priors
  target += normal_lpdf(b0 | 100, 5);
  target += normal_lpdf(b | 0, 5);

  target += exponential_lpdf(sigma | 0.5);
  target += normal_lpdf(rho_un | 0, 1);
  
  vector[N] eta;
  
  for(i in 1:Q){
    
    int N_id = e_id[i] - s_id[i] + 1;
    vector[N_id] mu = b0 + X[s_id[i]:e_id[i]] * b;
    vector[N_id] epsi;
    
    epsi[1] = y[s_id[i]] - mu[1];
    
    for(j in 2:N_id){
      epsi[j] = (y[s_id[i] + j - 1] - mu[j]);  
      mu[j] = mu[j] + rho*epsi[j-1];
    }
    
    eta[s_id[i]:e_id[i]] = mu;
    
  }
  
  target += normal_lpdf(y | eta, sigma);
  
  
}

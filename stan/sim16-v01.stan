data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;

  // pt level auc - assumes no missing in the 10 days...
  // actually, assumes that missing equates to zero, which would be wrong, clearly.
  vector<lower=0>[N] y;
  
  vector[2] pri_b_0;
  vector[2] pri_b;
  real pri_s;
  
  int prior_only;
}
parameters {
  real b_0;
  vector[P] b;
  real<lower=0> s;
}
model {

  target += normal_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  target += normal_lpdf(b | pri_b[1], pri_b[2]);
  target += exponential_lpdf(s | pri_s);

  if(!prior_only){
    target += normal_lpdf(y | b_0 + X*b, s);
  }
}
generated quantities{
  
}



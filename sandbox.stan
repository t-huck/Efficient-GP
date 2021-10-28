data {
  int<lower=0> N;
  int<lower=0> D;
  vector[D] x[N];
  vector[D] ell;
}


transformed data{
  vector[D] s[N];
  matrix[N,N] K;
  for (n in 1:N) {
    s[n] = x[n].* ell;
  }
  
  K = cov_exp_quad(s, 2.0, 1.0);
}

parameters {
  vector[N] mu;
  vector[N] y;
}


model {
  y ~ multi_normal(mu, K);
}


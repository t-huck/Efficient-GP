data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> M;
  real y[N];
  vector[D] x[N];
  vector[D] x_pred[M];

  vector[D] ell;
  real<lower=0> sf;
  real<lower=0> sn;
}

transformed data{
  vector[D] s[N];
  matrix[N,N] K;
  vector[D] s_pred[N];
  matrix[N, N] Kn;
  matrix[N, N] L;
  real alpha[N];
  matrix[N, M] Kmn;
  for (n in 1:N) {
    s[n] = x[n].* ell;
  }
  for (n in 1:N) {
    s_pred[n] = x_pred[n].* ell;
  }
  Kn =   cov_exp_quad(s, sf, 1.0)
                     + diag_matrix(rep_vector(sn, N));
  L = cholesky_decompose(Kn);
  alpha = algebra_solver(solver, 1, L', algebra_solver(solver, 1, L , y))
  Kmn =   cov_exp_quad(s[N]./ell, s_pred[N]./ell, sf, 1);
  matrix[M, M] Km =   cov_exp_quad(s_pred[N]./ell, sf, 1);
  matrix[N, M] V = algebra_solver(solver, 1, L, Kmn);
  real mu[N] = Kmn*alpha;
  matrix[M, M] S = Km - V'*V;
}

parameters {
}

model {
} 


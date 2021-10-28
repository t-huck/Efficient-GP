data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> M;
  vector[D] x[N];
  vector[D] xpred[M];

  vector[D] ell;
  real<lower=0> sf;
  real<lower=0> sn;
}

transformed data{
  matrix[N, N] Kn =   cov_exp_quad(x[N]./ell, sf, 1)
                     + diag_matrix(rep_vector(sn, N));
  matrix[N, N] L = cholesky_decompose(cov);
  real alpha[N] = algebra_solver(solver, 1, L', algebra_solver(solver, 1, L, y))
  matrix[N, Np] Kmn =   cov_exp_quad(x[N]./ell, xpred[N]./ell, sf, 1);
  matrix[Np, Np] Km =   cov_exp_quad(xpred[N]./ell, sf, 1);
  matrix[N, Np] V = algebra_solver(solver, 1, L, Kmn);
  real mu[N] = Kmn*alpha;
  matrix[Np, Np] S = Km - V'*V;
}

parameters {
}

model {
} 


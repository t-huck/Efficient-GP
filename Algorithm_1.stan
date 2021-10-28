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
  vector[N] alpha;
  matrix[N, M] Kmn;
  vector[N] y_vec;
  matrix[M, M] Km;
  matrix[N, M] V;
  vector[N] mu;
  matrix[M, M] S;
  matrix[M, M] Lm;
  for (n in 1:N) {
    s[n] = x[n].* ell;
    y_vec[n] = y[n];
  }
  for (n in 1:M) {
    s_pred[n] = x_pred[n].* ell;
  }
  Kn =   cov_exp_quad(s, sf, 1.0)
                     + diag_matrix(rep_vector(sn, N));
  L = cholesky_decompose(Kn);
  alpha = mdivide_right_tri_low(mdivide_left_tri_low(L, y_vec)', L)';
  Kmn =   cov_exp_quad(s, s_pred, sf, 1.0);
  Km =   cov_exp_quad(s_pred, sf, 1.0);
  V = mdivide_left_tri_low(L, Kmn);
  mu = Kmn*alpha;
  S = Km - V'*V;
  Lm = cholesky_decompose(S);
}

parameters {
  vector[M] eta;
}

model {
  eta ~ std_normal();
} 

generated quantities{
  vector[M] y_pred;
  y_pred = Lm*eta + mu;
}



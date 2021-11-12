functions {
  //Returns Phi and Q for Gaussian w/ Matern(7/2) used for approximation
  //Algorithm assumes the data is unidimensional
  matrix getPhi(real d, real lengthscale){
    real l = sqrt(7)/lengthscale;
    matrix[4, 4] A = [[(l^3*d^3+3*l^2*d^2+6*l*d+6),  3*(l^2*d^3+2*l*d+2*d), 3*(l*d^3+d^2), d^3],
     [-l^4*d^3, -3*(l^3*d^3-l^2*d^2-2*l*d-2), -3*(l^2*d^3-2*l*d^2-2*d), -(l*d^3-3*d^2)],
    [(l^6*d^3-6*l^5*d^2+6*l^4*d),3*(l^4*d^3-4*l^3*d^2),3*(l^3*d^3-5*l^2*d^2+2*l*d+2),(l^2*d^3-6*l*d^2+6*d)],  
    [-(l^6*d^3-6*l^5*d^2+6*l^4*d), -3*(l^5*d^3-7*l^4*d^2+8*l^3*d),  
    -3*(d^4*d^3-8*l^3*d^2+12*l^2*d), -(l^3*d^3-9*l^2*d^2+18*l*d-6)]];
    print("l=",l);
    print("d=",d);
    print("A=", A);
    return A*exp(-l*d)/6;
  }
  
  matrix getQ(real d, real lengthscale, real sf){
    real l = sqrt(7)/lengthscale;
    real q = sf*l^7*32/5;
    matrix[4,4] Q = [[-(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 + 12 * exp(-2 * l * d) * d ^ 5 * l ^ 5 + 30 * exp(-2 * l * d) * d ^ 4 * l ^ 4 + 60 * exp(-2 * l * d) * d ^ 3 * l ^ 3 + 90 * exp(-2 * l * d) * d ^ 2 * l ^ 2 + 90 * exp(-2 * l * d) * d * l + 45 * exp(-2 * l * d) - 45) / l ^ 7 / 288,
    d ^ 6 * exp(-2 * l * d) / 72,
    -(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 12 * exp(-2 * l * d) * d ^ 5 * l ^ 5 - 6 * exp(-2 * l * d) * d ^ 4 * l ^ 4 - 12 * exp(-2 * l * d) * d ^ 3 * l ^ 3 - 18 * exp(-2 * l * d) * d ^ 2 * l ^ 2 - 18 * exp(-2 * l * d) * d * l - 9 * exp(-2 * l * d) + 9) / l ^ 5 / 288,
    exp(-2 * l * d) * d ^ 4 * (d ^ 2 * l ^ 2 - 6 * l * d + 3) / 72],
    [d ^ 6 * exp(-2 * l * d) / 72,
    -(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 12 * exp(-2 * l * d) * d ^ 5 * l ^ 5 + 6 * exp(-2 * l * d) * d ^ 4 * l ^ 4 + 12 * exp(-2 * l * d) * d ^ 3 * l ^ 3 + 18 * exp(-2 * l * d) * d ^ 2 * l ^ 2 + 18 * exp(-2 * l * d) * d * l + 9 * exp(-2 * l * d) - 9) / l ^ 5 / 288,
    d ^ 4 * exp(-2 * l * d) * (d ^ 2 * l ^ 2 - 6 * l * d + 9) / 72,
    -(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 36 * exp(-2 * l * d) * d ^ 5 * l ^ 5 + 90 * exp(-2 * l * d) * d ^ 4 * l ^ 4 - 60 * exp(-2 * l * d) * d ^ 3 * l ^ 3 - 18 * exp(-2 * l * d) * d ^ 2 * l ^ 2 - 18 * exp(-2 * l * d) * d * l - 9 * exp(-2 * l * d) + 9) / l ^ 3 / 288],
    [-(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 12 * exp(-2 * l * d) * d ^ 5 * l ^ 5 - 6 * exp(-2 * l * d) * d ^ 4 * l ^ 4 - 12 * exp(-2 * l * d) * d ^ 3 * l ^ 3 - 18 * exp(-2 * l * d) * d ^ 2 * l ^ 2 - 18 * exp(-2 * l * d) * d * l - 9 * exp(-2 * l * d) + 9) / l ^ 5 / 288,
    d ^ 4 * exp(-2 * l * d) * (d ^ 2 * l ^ 2 - 6 * l * d + 9) / 72,
    -(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 36 * exp(-2 * l * d) * d ^ 5 * l ^ 5 + 102 * exp(-2 * l * d) * d ^ 4 * l ^ 4 - 84 * exp(-2 * l * d) * d ^ 3 * l ^ 3 + 18 * exp(-2 * l * d) * d ^ 2 * l ^ 2 + 18 * exp(-2 * l * d) * d * l + 9 * exp(-2 * l * d) - 9) / l ^ 3 / 288,
    d ^ 2 * exp(-2 * l * d) * (d ^ 4 * l ^ 4 - 12 * d ^ 3 * l ^ 3 + 48 * d ^ 2 * l ^ 2 - 72 * l * d + 36) / 72],
    [exp(-2 * l * d) * d ^ 4 * (d ^ 2 * l ^ 2 - 6 * l * d + 3) / 72,
     -(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 36 * exp(-2 * l * d) * d ^ 5 * l ^ 5 + 90 * exp(-2 * l * d) * d ^ 4 * l ^ 4 - 60 * exp(-2 * l * d) * d ^ 3 * l ^ 3 - 18 * exp(-2 * l * d) * d ^ 2 * l ^ 2 - 18 * exp(-2 * l * d) * d * l - 9 * exp(-2 * l * d) + 9) / l ^ 3 / 288,
    d ^ 2 * exp(-2 * l * d) * (d ^ 4 * l ^ 4 - 12 * d ^ 3 * l ^ 3 + 48 * d ^ 2 * l ^ 2 - 72 * l * d + 36) / 72,
    -(4 * exp(-2 * l * d) * d ^ 6 * l ^ 6 - 60 * exp(-2 * l * d) * d ^ 5 * l ^ 5 + 318 * exp(-2 * l * d) * d ^ 4 * l ^ 4 - 708 * exp(-2 * l * d) * d ^ 3 * l ^ 3 + 666 * exp(-2 * l * d) * d ^ 2 * l ^ 2 - 198 * exp(-2 * l * d) * d * l + 45 * exp(-2 * l * d) - 45) /(l*288)]];
    return q*Q;
  }
  
  matrix getV(real x, real lengthscale, real sf){
    real l = sqrt(7)/lengthscale;
    real q = sf*l^7*32/5;
    matrix[4, 4] V = [[-(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 + 12 * exp(-2 * l * x) * l ^ 5 * x ^ 5 + 30 * exp(-2 * l * x) * l ^ 4 * x ^ 4 + 60 * exp(-2 * l * x) * l ^ 3 * x ^ 3 + 90 * exp(-2 * l * x) * l ^ 2 * x ^ 2 + 90 * exp(-2 * l * x) * l * x + 45 * exp(-2 * l * x) - 45) / l ^ 7 / 288,
    exp(-2 * l * x) * x ^ 6 / 72,
    -(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 12 * exp(-2 * l * x) * l ^ 5 * x ^ 5 - 6 * exp(-2 * l * x) * l ^ 4 * x ^ 4 - 12 * exp(-2 * l * x) * l ^ 3 * x ^ 3 - 18 * exp(-2 * l * x) * l ^ 2 * x ^ 2 - 18 * exp(-2 * l * x) * l * x - 9 * exp(-2 * l * x) + 9) / l ^ 5 / 288,
    x ^ 4 * (l ^ 2 * x ^ 2 - 6 * l * x + 3) * exp(-2 * l * x) / 72],
    [exp(-2 * l * x) * x ^ 6 / 72,
    -(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 12 * exp(-2 * l * x) * l ^ 5 * x ^ 5 + 6 * exp(-2 * l * x) * l ^ 4 * x ^ 4 + 12 * exp(-2 * l * x) * l ^ 3 * x ^ 3 + 18 * exp(-2 * l * x) * l ^ 2 * x ^ 2 + 18 * exp(-2 * l * x) * l * x + 9 * exp(-2 * l * x) - 9) / l ^ 5 / 288,
    (l * x - 3) ^ 2 * exp(-2 * l * x) * x ^ 4 / 72,
    -(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 36 * exp(-2 * l * x) * l ^ 5 * x ^ 5 + 90 * exp(-2 * l * x) * l ^ 4 * x ^ 4 - 60 * exp(-2 * l * x) * l ^ 3 * x ^ 3 - 18 * exp(-2 * l * x) * l ^ 2 * x ^ 2 - 18 * exp(-2 * l * x) * l * x - 9 * exp(-2 * l * x) + 9) / l ^ 3 / 288],
    [-(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 12 * exp(-2 * l * x) * l ^ 5 * x ^ 5 - 6 * exp(-2 * l * x) * l ^ 4 * x ^ 4 - 12 * exp(-2 * l * x) * l ^ 3 * x ^ 3 - 18 * exp(-2 * l * x) * l ^ 2 * x ^ 2 - 18 * exp(-2 * l * x) * l * x - 9 * exp(-2 * l * x) + 9) / l ^ 5 / 288,
    (l * x - 3) ^ 2 * exp(-2 * l * x) * x ^ 4 / 72,
    -(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 36 * exp(-2 * l * x) * l ^ 5 * x ^ 5 + 102 * exp(-2 * l * x) * l ^ 4 * x ^ 4 - 84 * exp(-2 * l * x) * l ^ 3 * x ^ 3 + 18 * exp(-2 * l * x) * l ^ 2 * x ^ 2 + 18 * exp(-2 * l * x) * l * x + 9 * exp(-2 * l * x) - 9) / l ^ 3 / 288,
    (l ^ 2 * x ^ 2 - 6 * l * x + 6) ^ 2 * exp(-2 * l * x) * x ^ 2 / 72],
    [x ^ 4 * (l ^ 2 * x ^ 2 - 6 * l * x + 3) * exp(-2 * l * x) / 72,
    -(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 36 * exp(-2 * l * x) * l ^ 5 * x ^ 5 + 90 * exp(-2 * l * x) * l ^ 4 * x ^ 4 - 60 * exp(-2 * l * x) * l ^ 3 * x ^ 3 - 18 * exp(-2 * l * x) * l ^ 2 * x ^ 2 - 18 * exp(-2 * l * x) * l * x - 9 * exp(-2 * l * x) + 9) / l ^ 3 / 288,
    (l ^ 2 * x ^ 2 - 6 * l * x + 6) ^ 2 * exp(-2 * l * x) * x ^ 2 / 72,
    -(4 * exp(-2 * l * x) * l ^ 6 * x ^ 6 - 60 * exp(-2 * l * x) * l ^ 5 * x ^ 5 + 318 * exp(-2 * l * x) * l ^ 4 * x ^ 4 - 708 * exp(-2 * l * x) * l ^ 3 * x ^ 3 + 666 * exp(-2 * l * x) * l ^ 2 * x ^ 2 - 198 * exp(-2 * l * x) * l * x + 45 * exp(-2 * l * x) - 45) / (l* 288)]];
    return q*V;
  }
}

data {
  int<lower=1> N;
  vector[N] x_raw;
  vector[N] y_raw;
  vector[N] train;
  real ell;
  real<lower=0> sf;
  real<lower=0> sn;
}

transformed data{
  vector[N] x = (x_raw - mean(x_raw)/sd(x_raw));
  vector[N] y = (y_raw - mean(y_raw)/sd(y_raw));
}

parameters {
  
}

model {
  matrix[4,4] Phi_mat[N]; //Ob1
  matrix[4,4] Q_mat[N]; //Ob1
  matrix[4,4] P_mat[N]; 
  vector[4] G_mat[N]; // Indexing **not** offset by 1 from original paper
  matrix[4,4] V[N+1]; //Ob1
  matrix[4,4] L[N];
  vector[4] mu[N+1]; //Offset by 1 
  real v_star[N];
  real mu_star[N];
  matrix[4, 4] Ezz_same[N];
  matrix[4, 4] Ezz_step[N];
  row_vector[4] H = [1,0,0,0];
  mu[1] = [0,0,0,0]';
  V[1] = getV(x[1], ell, sf);
  print("V[1]=", V[1]);
  for(n in 2:N+1){
    if(n==2){
      Phi_mat[n-1] = getPhi(50, ell);
      Q_mat[n-1] =  getQ(50, ell, sf);
    }
    else{
      Phi_mat[n-1] = getPhi(x[n-1]-x[n-2], ell);
      Q_mat[n-1] =  getQ(x[n-1]-x[n-2], ell, sf);
    }
    print("Q_mat[",n-1,"]=", Q_mat[n-1]);
    print("Phi_mat[", n-1, "]=", Phi_mat[n-1]);
    print("V[",n-1,"]=",V[n-1]);
    print("mu[",n-1,"]=", mu[n-1]);
    P_mat[n-1] = Phi_mat[n-1]*V[n-1]*Phi_mat[n-1]'+Q_mat[n-1];
    print("P_mat[n-1]=", P_mat[n-1]);
    G_mat[n] = P_mat[n-1]*H'/(H*P_mat[n-1]*H' + sn);
    print("G_mat[",n-1,"]=", G_mat[n-1]);
    if(train[n-1]){
      target += normal_lpdf(y[n-1]| H*Phi_mat[n-1]*mu[n-1], H*P_mat[n-1]*H'+sn);
      mu[n] = Phi_mat[n-1]*mu[n-1] + G_mat[n]*(y[n]-H*Phi_mat[n-1]*mu[n-1]);
      V[n] = P_mat[n-1] - G_mat[n]*H*P_mat[n-1];
    }
    else{
        mu[n] = Phi_mat[n-1]*mu[n-1];
        V[n] = P_mat[n-1];
    }
  }
    mu_star[N]=H*mu[N+1];
    v_star[N] = H*V[N+1]*H';
    Ezz_same[N] = V[N+1] + mu[N+1]*mu[N+1]';
    Ezz_step[N] = (diag_matrix(rep_vector(1.0, 4)) -G_mat[N]*H)*Phi_mat[N]*V[N+1];
    L[N] = V[N+1]*Phi_mat[N]*inverse(P_mat[N]);
    for(t in N-1:1){
      L[t] = V[t+1]*Phi_mat[t]*inverse(P_mat[t]);
      mu[t] += L[t]*(mu[t+1]-Phi_mat[t]*mu[t]);
      mu_star[t] = H*mu[t];
      V[t] += L[t]*(V[t+1]-P_mat[t])*L[t]';
      v_star[t] = H*V[t]*H';
      Ezz_same[t] = V[t] + mu[t]*mu[t]';
      Ezz_step[t] = V[t+1]*L[t]' + L[t+1]*(Ezz_step[t-1] - Phi_mat[t+1]*V[t+1])*L[t]';
    }
    
} 

generated quantities{

}



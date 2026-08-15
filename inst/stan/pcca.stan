data {
  int<lower=3> n;
  int<lower=1> p;
  int<lower=1> q;
  int<lower=1> k;
  matrix[n, p] X;
  matrix[n, q] Y;
}
parameters {
  matrix[p, k] Wx;
  matrix[q, k] Wy;
  matrix[n, k] z;
  vector<lower=0>[p] sigma_x;
  vector<lower=0>[q] sigma_y;
}
model {
  to_vector(Wx) ~ normal(0, 1);
  to_vector(Wy) ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  sigma_x ~ normal(0, 1);
  sigma_y ~ normal(0, 1);
  for (i in 1:n) {
    to_vector(X[i]) ~ normal(Wx * to_vector(z[i]), sigma_x);
    to_vector(Y[i]) ~ normal(Wy * to_vector(z[i]), sigma_y);
  }
}

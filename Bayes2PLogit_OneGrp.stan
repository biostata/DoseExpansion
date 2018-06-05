data {
  int<lower=0> n;
  int<lower=0,upper=1> y[n];
  real<lower=1> x[n];
}
parameters {
  real mu;
  real<lower=0> beta;
}
model {
  vector[n] x_beta;
  mu ~ cauchy(0,10);
  beta ~ cauchy(0,2.5);
  for(i in 1:n) 
    x_beta[i] = mu + x[i] * beta;
  y ~ bernoulli_logit(x_beta);
}
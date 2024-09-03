data {
  int<lower=0> N;     // size of data
  int<lower=0> steps; // number of steps per track
  vector[N] t;        // time step
  vector[N] theta;    // step direction
  vector[N] theta1;   // previous step direction
  vector[N] gamma;    // bearing to target (not used in process1)
  vector[N] dist;     // distance to target (not used in process1)
}
parameters {
  real<lower=0,upper=steps> t_hat; // time step at prey detection
  real<lower=0> kappa1;        // kappa prior to t_hat
  real<lower=0> dkappa;        // increase in kappa after t_hat
}
model {
  real beta = 0;
  vector[N] mu_x, mu_y;
  real mu;
  // Likelihood of the walk given the parameters
  mu_x = cos(theta1) * (1 - beta) + cos(gamma) * beta;
  mu_y = sin(theta1) * (1 - beta) + sin(gamma) * beta;
  for (i in 1:N) {
    // atan2 is vectorized in stan >= 2.34, but rstan is currently on 2.32.2
    mu = atan2(mu_y[i], mu_x[i]);
    if (t[i] < t_hat) {
      theta[i] ~ von_mises(mu, kappa1);
    } else {
      theta[i] ~ von_mises(mu, kappa1 + dkappa);
    }
  }
  // Priors
  t_hat ~ uniform(1, steps);
  kappa1 ~ gamma(2, 0.05);
  dkappa ~ gamma(2, 0.1);
}

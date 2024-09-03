library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

theme_set(theme_classic())
set.seed(1001)

# Starting parameters
x0 <- y0 <- b0 <- 0
nsteps <- 1000
sigma <- 1
beta <- 0.1
# Attractor
xn <- yn <- nsteps / 2
attractor <- tibble(x = xn, y = yn)

# Continuous time correlated random walk ----------------------------------

# v_c(t+1) = gamma_c + exp(-beta)(v_c(t) - gamma_c) + epsilon_c(1); var(epsilon)=sigma^2(1-exp(-2beta))/(2beta)

vx <- vy <- numeric(nsteps)
ctcrw_sd <- sqrt(sigma^2 * (1 - exp(-2 * beta)) / (2 * beta))
epsilon_x <- rnorm(nsteps, mean = 0, sd = ctcrw_sd)
epsilon_y <- rnorm(nsteps, mean = 0, sd = ctcrw_sd)
vx[1] <- epsilon_x[1]
vy[1] <- epsilon_y[1]
for (t in 2:nsteps) {
  vx[t] <- exp(-beta) * vx[t - 1] + epsilon_x[t]
  vy[t] <- exp(-beta) * vy[t - 1] + epsilon_y[t]
}

x <- x0 + cumsum(vx)
y <- y0 + cumsum(vy)

crw <- tibble(
  t = seq(nsteps),
  vx = vx,
  vy = vy,
  x = x,
  y = y
)

ggplot(crw, aes(x, y)) +
  geom_path() +
  geom_point(aes(color = factor(t)), slice(crw, 1, nsteps), size = 2) +
  geom_point(data = attractor, size = 2, color = "orange") +
  coord_fixed() +
  theme(legend.position = "none")


# Continuous time biased random walk --------------------------------------

vx <- vy <- x <- y <- numeric(nsteps)
ctcrw_sd <- sqrt(sigma^2 * (1 - exp(-2 * beta)) / (2 * beta))
epsilon_x <- rnorm(nsteps, mean = 0, sd = ctcrw_sd)
epsilon_y <- rnorm(nsteps, mean = 0, sd = ctcrw_sd)
vx[1] <- epsilon_x[1]
vy[1] <- epsilon_y[1]
x[1] <- x0
y[1] <- y0
for (t in 2:nsteps) {
  bearing <- atan2(yn - y[t - 1], xn - x[t - 1])
  gamma_x <- cos(bearing)
  gamma_y <- sin(bearing)
  vx[t] <- gamma_x + exp(-beta) * (vx[t - 1] - gamma_x) + epsilon_x[t]
  vy[t] <- gamma_y + exp(-beta) * (vy[t - 1] - gamma_y) + epsilon_y[t]
  x[t] <- x[t - 1] + vx[t]
  y[t] <- y[t - 1] + vy[t]
}

brw <- tibble(
  t = seq(nsteps),
  vx = vx,
  vy = vy,
  x = x,
  y = y
)

ggplot(brw, aes(x, y)) +
  geom_path() +
  geom_point(aes(color = factor(t)), slice(brw, 1, nsteps), size = 2) +
  geom_point(data = attractor, size = 2, color = "orange") +
  coord_fixed() +
  theme(legend.position = "none")

# Continuous time biased random walk 2 ------------------------------------

vx <- vy <- x <- y <- numeric(nsteps)
ctcrw_sd <- sqrt(sigma^2 * (1 - exp(-2 * beta)) / (2 * beta))
epsilon_x <- rnorm(nsteps, mean = 0, sd = ctcrw_sd)
epsilon_y <- rnorm(nsteps, mean = 0, sd = ctcrw_sd)
vx[1] <- epsilon_x[1]
vy[1] <- epsilon_y[1]
x[1] <- x0
y[1] <- y0
for (t in 2:nsteps) {
  bearing <- atan2(yn - y[t - 1], xn - x[t - 1])
  distance <- sqrt((yn - y[t - 1])^2 + (xn - x[t - 1])^2)
  gamma_x <- cos(bearing)
  gamma_y <- sin(bearing)
  vx[t] <- gamma_x + exp(-beta) * (vx[t - 1] - gamma_x) + epsilon_x[t]
  vy[t] <- gamma_y + exp(-beta) * (vy[t - 1] - gamma_y) + epsilon_y[t]
  x[t] <- x[t - 1] + vx[t]
  y[t] <- y[t - 1] + vy[t]
}

brw <- tibble(
  t = seq(nsteps),
  vx = vx,
  vy = vy,
  x = x,
  y = y
)

ggplot(brw, aes(x, y)) +
  geom_path() +
  geom_point(aes(color = factor(t)), slice(brw, 1, nsteps), size = 2) +
  geom_point(data = attractor, size = 2, color = "orange") +
  coord_fixed() +
  theme(legend.position = "none")

# Ok but what if the *acceleration* is what's correlated/attracted?
# Treat it like a force?

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

theme_set(theme_classic())
set.seed(1001)

# Forces involved: inertia (CRW), attraction (BRW), and drag (keep speeds sane)

# Number of steps
nsteps <- 100
# Location of attractor (north of origin)
xa <- 0
ya <- nsteps

# Utility functions
dot <- function(A, B) sum(A * B)
norm <- function(A) sqrt(dot(A, A))
turn_angle <- function(A, B) {
  acos(dot(A, B) / (norm(A) * norm(B)))
}

# CRW ---------------------------------------------------------------------

# Parameters
sigma <- 1    # Random walk variance^2
cd <- 0.01       # Drag coefficient

# Acceleration, velocity, location
ax <- ay <- vx <- vy <- x <- y <- numeric(nsteps)
# Derived values: turning angle (theta), bearing (phi), drag
theta <- phi <- dragx <- dragy <- numeric(nsteps)

# Initial conditions (accelerating 1 unit north, velocity and position are 0)
ax[1] <- 0
ay[1] <- 1
vx[1] <- vy[1] <- x[1] <- y[1] <- 0
theta[1] <- phi[1] <- dragx[1] <- dragy[1] <- NA

# Random walk errors
dax <- day <- rnorm(nsteps, 0, sigma)

# Run simulation
for (t in 2:nsteps) {
  # Add random walk
  axt <- ax[t - 1] + dax[t]
  ayt <- ay[t - 1] + day[t]
  # Add drag
  vt <- c(vx[t - 1], vy[t - 1])  # Velocity vector
  dragt <- -cd * vt * norm(vt) # Fd = cd * v * ||v||
  dragx[t] <- dragt[1]
  dragy[t] <- dragt[2]
  ax[t] <- axt + dragx[t]
  ay[t] <- ayt + dragy[t]
  # Velocity and location (then turning angle and bearing)
  vx[t] <- vx[t - 1] + ax[t]
  vy[t] <- vy[t - 1] + ay[t]
  x[t] <- x[t - 1] + vx[t]
  y[t] <- y[t - 1] + vy[t]
  if (t == 2) {
    theta[t] <- NA
  } else {
    A <- c(x[t - 1] - x[t - 2], y[t - 1] - y[t - 2])
    B <- c(x[t] - x[t - 1], y[t] - y[t - 1])
    theta[t] <- turn_angle(A, B)
  }
  phi[t] <- atan2(y[t] - y[t - 1], x[t] - x[t - 1])
}

crw <- tibble(t = seq(nsteps),
              ax, ay, vx, vy, x, y,
              theta, phi, dragx, dragy) %>%
  mutate(a = sqrt(ax^2 + ay^2),
         v = sqrt(vx^2 + vy^2),
         drag = sqrt(dragx^2 + dragy^2))

ggplot(crw, aes(x, y)) +
  geom_path()

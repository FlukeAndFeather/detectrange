library(patchwork)
library(tidyverse)

theme_set(theme_classic())
set.seed(1001)

# Number of steps
nsteps <- 1000
# Location of attractor (north of origin)
xa <- 0
ya <- nsteps

# theta and phi are turning angle and bearing to attractor, respectively

normalize_angle <- \(theta) (theta + pi) %% (2 * pi) - pi

# CRW ---------------------------------------------------------------------

# Parameters
k_crw <- 20

# Initial conditions
x <- y <- vx <- vy <- theta <- phi <- numeric(nsteps)
x[1] <- y[1] <- vx <- vy <- theta[1] <- phi[1] <- 0

# Random walk
dtheta <- CircStats::rvm(nsteps, 0, k_crw)

for (t in 2:nsteps) {
  theta[t] <- normalize_angle(theta[t - 1] + dtheta[t])
  vx[t] <- cos(theta[t])
  vy[t] <- sin(theta[t])
  x[t] <- x[t - 1] + vx[t]
  y[t] <- y[t - 1] + vy[t]
}

crw <- tibble(t = seq(nsteps), x, y, vx, vy, theta)
endpoints <- crw[c(1, nsteps),]
endpoints$name <- factor(c("start", "end"), levels = c("start", "end"))

# Track plot
ggplot(crw, aes(x, y)) +
  geom_path() +
  geom_point(aes(color = name), endpoints, size = 2) +
  theme(legend.title = element_blank(),
        legend.position = c(0.01, 1),
        legend.justification = c(0, 1))

# Bearing and distance to final point
crw <- mutate(crw,
              # bearing to final point
              phi = atan2(y[nsteps] - lag(y), x[nsteps] - lag(x)),
              # diff btw track bearing (theta) and bearing to final point (phi)
              delta = normalize_angle(theta - phi),
              dist = sqrt((y[nsteps] - y)^2 + (x[nsteps] - x)^2))
ggplot(crw, aes(t, delta)) +
  geom_line()

# Detection distance by threshold angle
detect <- tibble(
  delta_thr = seq(min(abs(crw$delta), na.rm = TRUE), pi, length.out = 100),
  detect = map_dbl(
    delta_thr,
    \(d) {

      max(rev(crw$dist)[cummax(abs(rev(crw$delta))) <= d], na.rm = TRUE)
    }
  ))

ggplot(crw, aes(dist, delta)) +
  geom_path() +
  geom_ribbon(aes(x = detect, ymin = -delta_thr, ymax = delta_thr),
              detect,
              inherit.aes = FALSE,
              color = NA, fill = "firebrick", alpha = 0.5) +
  scale_y_continuous(breaks = seq(-pi, pi, by = pi / 4),
                     labels = seq(-180, 180, by = 45))

# Basically, distance vs delta is not a function because animals can get closer
# and farther through time. You can do time vs delta, but then detection
# distance stops making sense because for a given delta there can be an earlier
# point in time where the whale was *closer* and didn't detect the source.

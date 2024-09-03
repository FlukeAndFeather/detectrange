library(patchwork)
library(tidyverse)

theme_set(theme_classic())
set.seed(1001)

#' Normalize angle to [-pi, pi)
#'
#' @param theta Angle(s) in radians
#'
#' @return Angle(s), theta2, that meet two conditions. 1) theta2 = theta +
#'   2*k*pi, where k is an integer. 2) theta2 is an element of the range [-pi,
#'   pi).
#'
#' @examples
#' normalize_angle(c(-3*pi/2, -pi/2, 0, pi, 3*pi/2))
#' # pi/2, -pi/2, 0, -pi, -pi/2
normalize_angle <- \(theta) {
  (theta + pi) %% (2 * pi) - pi
}

#' Circular weighted mean
#'
#' The circular weighted mean is calculated using vectors.
#' 1) Create unit vectors for the angles in a
#' 2) Lengthen the unit vectors by the weights in w
#' 3) Sum all the vectors
#' 4) The angle of the resultant vector is the circular weighted mean
#'
#' @param a Angles in radians
#' @param w Weights
#'
#' @return A circular mean of a weighted by w
#' @export
#'
#' @examples
#' # Circular mean of two angles on either side of pi will be pi, not 0
#' circ_wgt_mean(c(pi - 0.05, -pi + 0.05), c(1, 1))
#' # pi
#' mean(c(pi - 0.05, -pi + 0.05))
#' # 0
circ_wgt_mean <- \(a, w) {
  # Using vector sums
  ax <- cos(a) * w
  ay <- sin(a) * w
  x <- sum(ax)
  y <- sum(ay)
  atan2(y, x)
}

#' Simulate a biased, correlated random walk
#'
#' @param nsteps number of steps
#' @param attractor location of a point attractor (source of bias)
#' @param rho dispersion of the step distribution
#' @param beta relative influence of correlated and bias random walks (0 = CRW,
#'   1 = BRW)
#' @param history number of previous steps influencing next step
#' @param history_wgt weight of previous steps. Must be length `history`. First
#'   element weights t-1, second t-2, etc.
#'
#' @return a data frame with columns:
#'   * t timestep
#'   * x x-coordinate
#'   * y y-coordinate
#'   * dx horizontal step length
#'   * dy vertical step length
#'   * theta step bearing
#'   * gamma bearing to attractor
#'   * mu turning angle expectation (0 if beta = 0, gamma if beta = 1)
#' @export
#'
#' @examples
#' sim_bcrw(1000, attractor = c(1000, 1000), rho = 0.8, beta = 0.1)
#' sim_bcrw(1000, attractor = c(1000, 1000), rho = 0.8, beta = 0.1, history = 4)
#' sim_bcrw(1000, attractor = c(1000, 1000), rho = 0.8, beta = 0.1, history = 4)
sim_bcrw <- function(nsteps, attractor, rho, beta, history = 1, history_wgt = rep(1, history)) {
  # Initial conditions
  x <- y <- dx <- dy <- theta <- gamma <- mu <- numeric(nsteps)
  x[1] <- y[1] <- dx <- dy <- theta[1] <- mu[1] <- 0
  gamma[1] <- atan2(attractor[2] - y[1], attractor[1] - x[1])

  # Simulate walk
  for (t in 2:nsteps) {
    # Weight the history
    history_idx <- max(1, t - history):(t - 1)
    theta_n <- circ_wgt_mean(theta[history_idx],
                             history_wgt[1:length(history_idx)])
    mu[t] <- circ_wgt_mean(c(theta_n, gamma[t - 1]), c(1 - beta, beta))
    theta[t] <- CircStats::rwrpcauchy(1, mu[t], rho)
    dx[t] <- cos(theta[t])
    dy[t] <- sin(theta[t])
    x[t] <- x[t - 1] + dx[t]
    y[t] <- y[t - 1] + dy[t]
    gamma[t] <- atan2(attractor[2] - y[t], attractor[1] - x[t])
  }

  tibble(t = seq(nsteps), x, y, dx, dy, theta, gamma, mu)
}

#' Simulates a biased correlated random walk from a known endpoint
#'
#' @param nsteps
#' @param attractor
#' @param rho
#' @param beta
#' @param history
#' @param history_wgt
#'
#' @return
#' @export
#'
#' @examples
sim_bcrw2 <- function(nsteps, attractor, rho, beta, history = 1, history_wgt = rep(1, history)) {
  # Initial conditions
  x <- y <- dx <- dy <- theta <- gamma <- mu <- numeric(nsteps)
  dx[1] <- dy[1] <- theta[1] <- 0

  # End at the attractor
  x[1] <- attractor[1]
  y[1] <- attractor[2]
  gamma[1] <- 0

  # Simulate walk
  for (t in 2:nsteps) {
    # Weight the history
    history_idx <- max(1, t - history):(t - 1)
    theta_n <- circ_wgt_mean(theta[history_idx],
                             history_wgt[1:length(history_idx)])
    mu[t] <- circ_wgt_mean(c(theta_n, gamma[t - 1]), c(1 - beta, beta))
    theta[t] <- normalize_angle(CircStats::rwrpcauchy(1, mu[t], rho))
    # Reverse direction
    dx[t] <- cos(theta[t])
    dy[t] <- sin(theta[t])
    x[t] <- x[t - 1] - dx[t]
    y[t] <- y[t - 1] - dy[t]
    gamma[t] <- atan2(attractor[2] - y[t], attractor[1] - x[t])
  }

  tibble(t = nsteps:1,
         x = x - attractor[1],
         y = y - attractor[2],
         dx, dy, theta, gamma, mu,
         bearing_error = normalize_angle(gamma - theta),
         dist = sqrt(x^2 + y^2)) %>%
    arrange(t)
}

# Find detection distance by bearing threshold
detection_distance <- function(rw) {
  find_t_thr <- function(t, bearing_error, thr) {
    bearing_error2 <- abs(rev(bearing_error))
    t2 <- rev(t)
    bearing_error_max <- cummax(bearing_error2)
    map_dbl(thr, \(.thr) min(t2[bearing_error_max <= .thr]))
  }
  tibble(
    bearing_error_thr = seq(0, pi / 2, length.out = 100),
    t = find_t_thr(rw$t, rw$bearing_error, bearing_error_thr),
    dist = map_dbl(t, \(.t) rw$dist[.t])
  )
}

plot_path <- function(rw) {
  endpoints <- rw[c(1, nrow(rw)), ] %>%
    mutate(name = c("start", "end"))
  track <- ggplot(rw, aes(x, y)) +
    geom_path() +
    geom_point(aes(color = name), endpoints) +
    coord_fixed() +
    theme(legend.position = "top",
          legend.title = element_blank())
  target <- ggplot(rw, aes(t, bearing_error)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line() +
    scale_y_continuous("Bearing error",
                       breaks = seq(-pi, pi, by = pi / 4),
                       labels = seq(-180, 180, by = 45)) +
    expand_limits(y = c(-pi, pi))
  detect_dist <- detection_distance(rw)
  detect <- ggplot(detect_dist, aes(bearing_error_thr, dist)) +
    geom_line() +
    scale_x_continuous("Bearing threshold",
                       breaks = seq(0, pi/2, by = pi / 12),
                       labels = seq(0, 90, by = 15))

  cowplot::plot_grid(
    track,
    cowplot::plot_grid(target, detect, nrow = 1),
    ncol = 1
  )
}

crw <- sim_bcrw2(200, c(1000, 1000), rho = 0.95, beta = 0)
plot_path(crw)

# What's the distribution of dist x bearing threshold

sim_tail <- function(epsilon_thr, kappa, beta, max_t = 1e3) {
  # Initial settings
  # x,y     = position
  # theta   = step direction
  # gamma   = bearing to target
  # epsilon = error (difference of theta, gamma)
  x <- y <- theta <- gamma <- epsilon <- numeric(max_t)
  x[1] <- y[1] <- epsilon[1] <- 0
  theta[1] <- gamma[1] <- pi / 4

  t <- 2
  while (abs(epsilon[t - 1]) <= epsilon_thr && t < max_t) {
    # Update location from previous step direction
    x[t] <- x[t - 1] - cos(theta[t - 1])
    y[t] <- y[t - 1] - sin(theta[t - 1])
    # Draw next step
    gamma[t] <- atan2(-y[t], -x[t])
    mu <- circ_wgt_mean(c(theta[t - 1], gamma[t]), c(1 - beta, beta))
    theta[t] <- normalize_angle(CircStats::rvm(1, mu, kappa))
    epsilon[t] <- normalize_angle(gamma[t] - theta[t])
    t <- t + 1
  }

  head(tibble(t = seq(max_t), x, y, theta, gamma, epsilon), t - 1) %>%
    mutate(dist = sqrt(x^2 + y^2))
}

plot_tail <- function(t) {
  angles <- pivot_longer(t, c(theta, gamma),
                         names_to = "angle", values_to = "radians") %>%
    mutate(xend = x + 0.75 * cos(radians),
           yend = y + 0.75 * sin(radians))
  epsilon_labels <- pretty(seq(min(t$epsilon) * 180/pi,
                               max(t$epsilon) * 180/pi,
                               length.out = 5))
  epsilon_breaks <- epsilon_labels * pi/180
  cowplot::plot_grid(
    ggplot(t, aes(x, y)) +
      geom_path() +
      geom_segment(aes(xend = xend, yend = yend, color = angle), angles) +
      annotate("point", 0, 0, color = "firebrick") +
      coord_fixed() +
      theme(legend.title = element_blank(),
            legend.position = "top"),
    ggplot(t, aes(t,  epsilon)) +
      geom_hline(yintercept = 0, color = "gray", linetype = "dotted") +
      geom_line() +
      scale_y_continuous("epsilon (deg)",
                         breaks = epsilon_breaks,
                         labels = epsilon_labels),
    ncol = 1,
    rel_heights = c(2, 1)
  )
}

crw <- sim_tail(epsilon_thr = 20 * pi/180, kappa = 40, beta = 0)
plot_tail(crw)

brw <- sim_tail(epsilon_thr = 20 * pi/180, kappa = 40, beta = 1)
plot_tail(brw)

cbrw <- sim_tail(epsilon_thr = 20 * pi/180, kappa = 40, beta = 0.5)
plot_tail(cbrw)

rnavdist <- function(n, epsilon_thr, kappa, beta, max_t = 1e3) {
  sapply(seq(n), \(.n) {
    t <- sim_tail(epsilon_thr, kappa, beta, max_t)
    max(t$dist)
  })
}

n_epsilon <- 40
n_iter <- 1000
kappa <- 10
betas <- c(0, 0.5, 1)
distances <- lapply(betas, \(.beta) {
  expand_grid(
    epsilon = seq(0, pi/4, length.out = n_epsilon),
    i = seq(n_iter)
  ) %>%
    group_by(epsilon) %>%
    mutate(navdist = rnavdist(n(), epsilon[1], kappa, beta = .beta)) %>%
    ungroup() %>%
    mutate(beta = .beta)
}) %>%
  list_rbind()
dist_plot <- function(d) {
  epsilon_labels <- pretty(seq(0, max(d$epsilon) * 180 / pi, length.out = 5))
  epsilon_breaks <- epsilon_labels * pi / 180
  d %>%
    group_by(beta, epsilon) %>%
    summarize(navdist_mean = mean(navdist),
              navdist_lwr = quantile(navdist, 0.1),
              navdist_upr = quantile(navdist, 0.9),
              .groups = "drop") %>%
    mutate(beta = factor(beta)) %>%
    ggplot(aes(epsilon)) +
    geom_ribbon(aes(ymin = navdist_lwr, ymax = navdist_upr, fill = beta),
                color = NA, alpha = 0.5) +
    geom_line(aes(y = navdist_mean, color = beta), linewidth = 1.5) +
    scale_x_continuous(expression(epsilon),
                       breaks = epsilon_breaks,
                       labels = epsilon_labels) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(y = "Distance")
}
cowplot::ggdraw(dist_plot(distances) + theme(legend.position = "none")) +
  cowplot::draw_plot(dist_plot(filter(distances, epsilon < pi / 6)) +
                       theme(legend.position = "inside",
                             legend.position.inside = c(0.01, 0.99),
                             legend.justification = c(0, 1)),
                     0.1, 0.4, 0.6, 0.6)




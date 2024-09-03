library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

theme_set(theme_classic())
set.seed(1001)
normalize_angle <- \(x_crw) (x_crw + pi) %% (2 * pi) - pi

# Starting parameters
x0 <- y0 <- b0 <- 0
nsteps <- 1000
# Attractor
xn <- yn <- nsteps

# Correlated random walk parameter kappa
k_crw <- 20

# CRW steps
a <- CircStats::rvm(n = nsteps, mean = 0, k = k_crw)
b_crw <- cumsum(a) %% (2 * pi)
x_crw <- x0 + cumsum(cos(b_crw))
y_crw <- y0 + cumsum(sin(b_crw))
# Bearing to attractor and difference from animal's bearing
bn <- atan2(yn - y_crw, xn - x_crw)
db <- bn - b_crw

# Visualize results
track_crw <- tibble(t = seq(length(x_crw)), x = x_crw, y = y_crw)
angles_crw <- tibble(t = seq(length(x_crw)),
                     a = normalize_angle(a),
                     db = normalize_angle(db))

track_plot_crw <- ggplot(track_crw, aes(x_crw, y_crw)) +
  geom_path() +
  annotate("point", track_crw$x[1], track_crw$y[1],
           color = "firebrick", size = 2) +
  annotate("point", last(track_crw$x), last(track_crw$y),
           color = "cornflowerblue", size = 2) +
  annotate("point", xn, yn,
           color = "orange", size = 2)

angles_plot_crw <- ggplot(angles_crw, aes(t, a * 180 / pi)) +
  geom_line() +
  geom_point(alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous("a",
                     breaks = seq(-180, 180, by = 45),
                     limits = c(-180, 180))

attractor_plot_crw <- ggplot(angles_crw, aes(t, db * 180 / pi)) +
  geom_line() +
  geom_point(alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous("b", breaks = seq(-180, 180, by = 45))

density_plots_crw <- angles_crw %>%
  pivot_longer(c(a, db), names_to = "angle", values_to = "value") %>%
  ggplot(aes(value * 180 / pi, color = angle)) +
  geom_density() +
  theme(legend.title = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1))

track_plot_crw + (angles_plot_crw / attractor_plot_crw / density_plots_crw)

# Biased random walk

k_brw <- 20
b_brw <- CircStats::rvm(n = nsteps, mean = 0, k = k_brw)
x_brw <- y_brw <- numeric(nsteps)
for (i in seq(nsteps)) {
  if (i == 1) {
    xi <- x0
    yi <- y0
  } else {
    xi <- x_brw[i - 1]
    yi <- y_brw[i - 1]
  }

  b <- atan2(yn - yi, xn - xi) + b_brw[i]

  x_brw[i] <- xi + cos(b)
  y_brw[i] <- yi + sin(b)
}
dx <- diff(c(x0, x_brw))
dy <- diff(c(y0, y_brw))
a_brw <- diff(c(b_brw[1], atan2(dy, dx)))

track_brw <- tibble(t = seq(length(x_brw)), x_brw, y_brw)
angles_brw <- tibble(t = seq(length(x_brw)),
                     a = normalize_angle(a_brw),
                     db = normalize_angle(b_brw))

track_plot_brw <- ggplot(track_brw, aes(x_brw, y_brw)) +
  geom_path() +
  annotate("point", track_brw$x_brw[1], track_brw$y_brw[1],
           color = "firebrick", size = 2) +
  annotate("point", last(track_brw$x_brw), last(track_brw$y_brw),
           color = "cornflowerblue", size = 2) +
  annotate("point", xn, yn,
           color = "orange", size = 2)

angles_plot_brw <- ggplot(angles_brw, aes(t, a * 180 / pi)) +
  geom_line() +
  geom_point(alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous("a",
                     breaks = seq(-180, 180, by = 45),
                     limits = c(-180, 180))

attractor_plot_brw <- ggplot(angles_brw, aes(t, db * 180 / pi)) +
  geom_line() +
  geom_point(alpha = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous("b",
                     breaks = seq(-180, 180, by = 45),
                     limits = c(-180, 180))

density_plots_brw <- angles_brw %>%
  pivot_longer(c(a, db), names_to = "angle", values_to = "value") %>%
  ggplot(aes(value * 180 / pi, color = angle)) +
  geom_density() +
  theme(legend.title = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1))

track_plot_brw + (angles_plot_brw / attractor_plot_brw / density_plots_brw)



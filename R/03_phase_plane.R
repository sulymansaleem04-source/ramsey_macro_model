
# ==============================================================================
# Phase 3: Phase Plane & Vector Field Visualization
# File: R/03_phase_plane.R
# ==============================================================================

library(ggplot2)
library(deSolve)

source("R/01_model_equations.R")
source("R/02_steady_state.R")

# 1. Steady State Coordinates
ss <- compute_steady_state(params)
k_star <- ss$k_star
c_star <- ss$c_star

k_max <- k_star * 2.0
c_max <- c_star * 2.0

# 2. Grid for Vector Field
grid <- expand.grid(
  k = seq(0.1, k_max, length.out = 25),
  c = seq(0.1, c_max, length.out = 25)
)

# 3. Compute Velocities dk and dc
grid$dk <- grid$k^params$alpha - grid$c - (params$n + params$g + params$delta) * grid$k
marginal_prod <- params$alpha * (grid$k^(params$alpha - 1))
hurdle        <- params$rho + params$delta + params$theta * params$g
grid$dc       <- (grid$c / params$theta) * (marginal_prod - hurdle)

# 4. Normalize Arrows
grid$len     <- sqrt(grid$dk^2 + grid$dc^2)
scale_factor <- 0.12
grid$dk_norm <- (grid$dk / grid$len) * scale_factor
grid$dc_norm <- (grid$dc / grid$len) * scale_factor

# 5. Nullcline Curves
k_seq <- seq(0.01, k_max, length.out = 300)
c_nullcline_k <- k_seq^params$alpha - (params$n + params$g + params$delta) * k_seq
df_nullcline  <- data.frame(k = k_seq, c = c_nullcline_k)
df_nullcline  <- df_nullcline[df_nullcline$c >= 0, ]

# 6. Plotting
phase_plot <- ggplot() +
  geom_segment(
    data = grid,
    aes(x = k, y = c, xend = k + dk_norm, yend = c + dc_norm),
    arrow = arrow(length = unit(0.14, "cm")),
    color = "grey55", alpha = 0.7
  ) +
  geom_line(
    data = df_nullcline,
    aes(x = k, y = c, color = "k_dot = 0"),
    linewidth = 1.2
  ) +
  geom_vline(
    aes(xintercept = k_star, color = "c_dot = 0"),
    linewidth = 1.2, linetype = "dashed"
  ) +
  geom_point(aes(x = k_star, y = c_star), color = "darkred", size = 4) +
  annotate(
    "text", x = k_star + 0.15, y = c_star + 0.08,
    label = sprintf("E* (k*=%.2f, c*=%.2f)", k_star, c_star),
    color = "darkred", fontface = "bold", hjust = 0
  ) +
  scale_color_manual(
    name = "Nullclines",
    values = c("k_dot = 0" = "#1f77b4", "c_dot = 0" = "#2ca02c")
  ) +
  labs(
    title = "Phase Plane: Ramsey-Cass-Koopmans Dynamic System",
    subtitle = "Directional field displaying saddle-path topology around steady state E*",
    x = "Capital per effective worker (k)",
    y = "Consumption per effective worker (c)"
  ) +
  coord_cartesian(xlim = c(0, k_max), ylim = c(0, c_max)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 7. Save and Display
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/01_phase_plane.png", plot = phase_plot, width = 8, height = 6, dpi = 300)
print(phase_plot)


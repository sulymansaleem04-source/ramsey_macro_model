# ==============================================================================
# Phase 4: Numerical Shooting Algorithm for Saddle-Path Trajectories
# File: R/04_shooting_solver.R
# ==============================================================================

library(deSolve)
library(ggplot2)

source("R/01_model_equations.R")
source("R/02_steady_state.R")

ss <- compute_steady_state(params)
k_star <- ss$k_star
c_star <- ss$c_star

solve_ramsey_shooting <- function(k0, p, T_horizon = 40, max_iter = 60, tol = 1e-6) {
  
  c_low  <- 1e-4
  c_high <- k0^p$alpha
  
  time_steps <- seq(0, T_horizon, by = 0.1)
  best_path  <- NULL
  
  cat(sprintf("Initiating shooting algorithm for k0 = %.3f (Target k* = %.3f)...\n", k0, k_star))
  
  for (iter in 1:max_iter) {
    c_guess <- (c_low + c_high) / 2
    init_state <- c(k = k0, c = c_guess)
    
    sim <- tryCatch({
      ode(
        y = init_state,
        times = time_steps,
        func = ramsey_ode_system,
        parms = p,
        method = "rk4"
      )
    }, error = function(e) NULL)
    
    if (is.null(sim)) {
      c_high <- c_guess
      next
    }
    
    sim_df <- as.data.frame(sim)
    best_path <- sim_df
    
    collapsed <- any(sim_df$k <= 0.01) || any(is.na(sim_df$k))
    exploded  <- any(sim_df$k > (k_star * 1.5))
    
    if (collapsed) {
      c_high <- c_guess
    } else if (exploded) {
      c_low <- c_guess
    } else {
      k_terminal <- tail(sim_df$k, 1)
      if (abs(k_terminal - k_star) < tol || (c_high - c_low) < tol) {
        cat(sprintf("Converged at iteration %d! Optimal c0 = %.5f\n", iter, c_guess))
        return(sim_df)
      }
      if (k_terminal < k_star) {
        c_high <- c_guess
      } else {
        c_low <- c_guess
      }
    }
  }
  
  cat(sprintf("Finished %d iterations. Approximate c0 = %.5f\n", max_iter, c_guess))
  return(best_path)
}

# 1. Run Shooting Method for Developing Economy (k0 = 1.0 < k*)
k_initial <- 1.0
trajectory <- solve_ramsey_shooting(k0 = k_initial, p = params)

# 2. Re-create base phase plane
source("R/03_phase_plane.R")

# 3. Overlay Saddle-Path Trajectory
saddle_plot <- phase_plot +
  geom_path(
    data = trajectory,
    aes(x = k, y = c),
    color = "darkblue",
    linewidth = 1.5
  ) +
  geom_point(
    aes(x = k_initial, y = trajectory$c[1]),
    color = "darkblue",
    size = 3.5
  ) +
  annotate(
    "text",
    x = k_initial + 0.1,
    y = trajectory$c[1] - 0.05,
    label = sprintf("Initial Jump: c(0)=%.3f", trajectory$c[1]),
    color = "darkblue",
    fontface = "bold",
    hjust = 0
  ) +
  labs(
    title = "Ramsey Model: Saddle-Path Trajectory via Shooting Method",
    subtitle = sprintf("Convergence from initial capital k0 = %.2f to steady state E*", k_initial)
  )

# 4. Save Figure
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/02_saddle_trajectory.png", plot = saddle_plot, width = 8, height = 6, dpi = 300)
print(saddle_plot)
cat("\nPhase 4 complete! Trajectory plot saved to output/figures/02_saddle_trajectory.png\n")








































































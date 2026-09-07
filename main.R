# ==============================================================================
# Master Execution Pipeline: Ramsey-Cass-Koopmans Dynamic Model
# File: main.R
# ==============================================================================

cat("========================================================\n")
cat(" Running Ramsey-Cass-Koopmans Macro Dynamic Pipeline\n")
cat("========================================================\n\n")

# 1. Load Core Model Modules
source("R/01_model_equations.R")
source("R/02_steady_state.R")

# 2. Compute Steady State & Stability Diagnostics
ss <- compute_steady_state(params)
J  <- compute_jacobian(params, ss)
eigen_vals <- eigen(J)$values

cat("1. Steady State Results:\n")
cat(sprintf("   - Capital (k*):     %.4f\n", ss$k_star))
cat(sprintf("   - Consumption (c*): %.4f\n", ss$c_star))
cat(sprintf("   - Output (y*):      %.4f\n", ss$k_star^params$alpha))
cat(sprintf("   - Stable Root:      %.4f\n", min(Re(eigen_vals))))
cat(sprintf("   - Explosive Root:   %.4f\n\n", max(Re(eigen_vals))))

# Save steady-state table to CSV
ss_table <- data.frame(
  Metric = c("Steady_State_k", "Steady_State_c", "Steady_State_y", "Stable_Eigenvalue", "Explosive_Eigenvalue"),
  Value  = c(ss$k_star, ss$c_star, ss$k_star^params$alpha, min(Re(eigen_vals)), max(Re(eigen_vals)))
)
write.csv(ss_table, "output/steady_state_summary.csv", row.names = FALSE)
cat("-> Saved: output/steady_state_summary.csv\n")

# 3. Generate Phase Plane & Vector Field
source("R/03_phase_plane.R")
cat("-> Saved: output/figures/01_phase_plane.png\n")

# 4. Solve Dynamic Saddle Path via Shooting Method
source("R/04_shooting_solver.R")
write.csv(trajectory, "output/saddle_path_trajectory.csv", row.names = FALSE)
cat("-> Saved: output/saddle_path_trajectory.csv\n")
cat("-> Saved: output/figures/02_saddle_trajectory.png\n\n")

cat("Pipeline completed successfully! All assets are saved in output/\n")




































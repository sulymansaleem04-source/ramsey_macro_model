# Calibrated Baseline Parameters
params <- list(
  alpha = 0.33,  # Capital share of output
  rho   = 0.04,  # Discount rate / impatience
  theta = 2.00,  # Risk aversion
  delta = 0.05,  # Depreciation rate
  n     = 0.01,  # Population growth
  g     = 0.02   # Technological progress rate
)

# Function 1: Compute Steady State Values
compute_steady_state <- function(p) {
  k_star <- (p$alpha / (p$rho + p$delta + p$theta * p$g))^(1 / (1 - p$alpha))
  c_star <- k_star^p$alpha - (p$n + p$g + p$delta) * k_star
  return(list(k_star = k_star, c_star = c_star))
}

# Function 2: Compute the Jacobian Matrix at Steady State
compute_jacobian <- function(p, ss) {
  J11 <- p$rho - p$n + (p$theta - 1) * p$g
  J12 <- -1
  J21 <- (ss$c_star * p$alpha * (p$alpha - 1) * (ss$k_star^(p$alpha - 2))) / p$theta
  J22 <- 0
  
  J <- matrix(c(J11, J21, J12, J22), nrow = 2, ncol = 2)
  rownames(J) <- c("k_dot", "c_dot")
  colnames(J) <- c("k", "c")
  return(J)
}























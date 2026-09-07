
library(deSolve)

ramsey_ode_system <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    k_val <- max(k, 1e-6)
    c_val <- max(c, 1e-6)
    
    y <- k_val^alpha
    dk <- y - c_val - (n + g + delta) * k_val
    
    marginal_product_k <- alpha * (k_val^(alpha - 1))
    hurdle_rate        <- rho + delta + theta * g
    dc <- (c_val / theta) * (marginal_product_k - hurdle_rate)
    
    return(list(c(dk, dc)))
  })
}


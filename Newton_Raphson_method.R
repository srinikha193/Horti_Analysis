px <- function(r, a, b, c, d, e, f) {
  a * ((4 - r^2) / (2 * r)) +
    b * (2 / r) +
    c * ((2 * r - 1) / (r^2 - r + 4)) +
    d * ((1 - 2 * r) / (r - r^2)) +
    e * ((-2 * (r - 1)) / (-r^2 + 2 * r + 3)) +
    f * ((2 * (r - 1)) / (1 - r)^2)
}
dpdx <- function(r, a, b, c, d, e, f) {
  a * ((8 - 2 * r^2 + 4) / (4 - r^2)^2) +
    b * (-2 / r^2) +
    c * ((-2 * r^2 + 2 * r + 7) / (r^2 - r + 4)^2) +
    d * ((-2 * r^2 + 2 * r - 1) / (r - r^2)^2) +
    e * ((-2 * r^2 + 4 * r - 10) / (-r^2 + 2 * r + 3)^2) +
    f * (-2 / (r^2 - 2 * r + 1))
}
# Set convergence tolerance and maximum iterations globally
eps <- 0.00005  # Convergence tolerance
max_iterations <- 50  # Maximum number of iterations

# Define the Newton-Raphson function to find the root
newton_raphson <- function(a, b, c, d, e, f, initial_guess = 0.5) {
  r <- initial_guess  # Start with an initial guess
  n <- 0  # Iteration counter
  
  # Print header for iteration details
  cat(sprintf("%-5s %-15s %-15s %-15s %-15s %-15s %s\n", "n", "r", "px(r)", "dpdx(r)", "r_new", "px(r_new)", "Update"))
  
  repeat {
    # Calculate dpdx and px at the current r value
    dpdx_val <- dpdx(r, a, b, c, d, e, f)
    px_val <- px(r, a, b, c, d, e, f)
    
    # Check for division by zero or invalid operations
    if (is.na(dpdx_val) || dpdx_val == 0 || is.nan(dpdx_val) || is.infinite(dpdx_val)) {
      cat("Derivative is zero or undefined. Stopping iteration.\n")
      return(NA)
    }
    
    # Newton's method update
    r_new <- r - (px_val / dpdx_val)
    px_new <- px(r_new, a, b, c, d, e, f)
    n <- n + 1
    
    # Print the iteration details
    cat(sprintf("%-5d %-15.10f %-15.10f %-15.10f %-15.10f %-15.10f %s\n",
                n, r, px_val, dpdx_val, r_new, px_new,
                ifelse(abs(px_new) < eps, "This is Root", "r = r_new")))
    
    # Check for convergence
    if (abs(px_new) < eps || n > max_iterations) {
      cat("Root is:", r_new, "\n")
      return(r_new)
    }
    
    # Update r for the next iteration
    r <- r_new
  }
}

# Example inputs for the coefficients
a <- maa  # Replace with your value
b <- mab  # Replace with your value
c <- mha  #
d <- mhb 
e <- mba  
f <- mbb  

# Run the Newton-Raphson method
root <- newton_raphson(a, b, c, d, e, f)
cat("Calculated Root:", root, "\n")

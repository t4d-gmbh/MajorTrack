
simulate_interactions <- function(N, probability_matrices, time_steps) {
  interaction_matrices <- list()
  
  step_idx <- 1
  for (m in seq_along(probability_matrices)) {
    for (t in 1:time_steps[m]) {
      interaction_matrix <- matrix(rbinom(N * N, 1, probability_matrices[[m]]), nrow=N, ncol=N)
      diag(interaction_matrix) <- 0  # No self-interactions
      interaction_matrices[[step_idx]] <- interaction_matrix
      step_idx <- step_idx + 1
    }
  }
  
  return(interaction_matrices)
}

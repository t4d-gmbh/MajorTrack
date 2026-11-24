# Minimal working example

# Define the problem parameters -----
N <- 10  # Number of agents
M <- 5   # Number of probability matrices
K <- 3   # Number of subgroups
subset_sizes <- c(3, 4, 3)  # Sizes of each subset
assignments_list <- lapply(1:M, function(x) generate_partitions(N, K, subset_sizes))

# Define the operations --------
operations <- list(
  list(type="merge", groups=c(1, 2)),
  list(type="split", group=3, sizes=c(2, 1)),
  list(type="drift", groups=c(1, 3), size=1)
)

# Build a starting assignment and modify it according to the operations -----
for (m in 2:M) {
  assignments_list[[m]] <- modify_partitions(assignments_list[[m - 1]], operations[[m %% length(operations) + 1]])
}

# Build the probability matrices ----------
p_intra_list <- replicate(M, c(0.8, 0.7, 0.9), simplify=FALSE)  # Different intra-group probabilities for each matrix
p_inter_list <- replicate(M, matrix(c(
  0.0, 0.2, 0.1,
  0.2, 0.0, 0.15,
  0.1, 0.15, 0.0
), nrow=K, byrow=TRUE), simplify=FALSE)  # Different inter-group probabilities

matrices <- generate_probability_matrices(N, M, assignments_list, p_intra_list, p_inter_list)

# Perform the simulations --------------
time_steps <- c(3, 2, 4, 3, 5)  # Number of time steps for each probability matrix
interaction_matrices <- simulate_interactions(N, matrices, time_steps)
print(interaction_matrices)


generate_probability_matrices <- function(N, M, assignments_list, p_intra_list, p_inter_list) {
  matrices <- list()
  #browser()
  for (m in 1:M) {
    assignments <- assignments_list[[m]]
    p_intra <- p_intra_list[[m]]
    p_inter <- p_inter_list[[m]]
    
    q_ab <- matrix(0, nrow=N, ncol=N)
    
    for (i in 1:N) {
      for (j in 1:N) {
        if (assignments[i] == assignments[j]) {
          q_ab[i, j] <- p_intra[assignments[i]]
        } else {
          q_ab[i, j] <- p_inter[assignments[i], assignments[j]]
        }
      }
    }
    
    matrices[[m]] <- q_ab
  }
  
  return(matrices)
}

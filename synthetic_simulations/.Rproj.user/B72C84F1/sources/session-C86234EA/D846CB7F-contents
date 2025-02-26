generate_partitions <- function(N, K, subset_sizes, assignments=NULL) {
  # Given a set of N agents and K groups of sizes specified in a vector
  # of integers such that sum(subset_sizes) = N, it randomly generate
  # a vector associating to each agent a group membership id in the
  # vector assignments.
  if (sum(subset_sizes) != N) {
    stop("Sum of subset sizes must equal N.")
  }
  
  if (is.null(assignments)) {
    assignments <- unlist(lapply(1:K, function(k) rep(k, subset_sizes[k])))
    assignments <- sample(assignments)
  }
  
  return(assignments)
}

modify_partitions <- function(assignments, operation) {
  # Modifies a vector of assignment according to rules specified in a list 
  # termed `operation`. The operation list should be specified as follows:
  #  * Merging groups: It merges the groups specified into a single group. The list
  #    should contain:
  #      * `operation$type = "merge"`. A string specifying the type of operation
  #      * `operation$groups`: A vector of integers specifing the groups to be merged into the
  #         one specified as first element of the vector.
  #  * Splitting a group: It split a groups into two or more subgroups of specified size.      should contain:
  #         * `operation$type = "split"`. A string specifying the type of operation
  #         * `operation$group`: An integer specifying the group to split.
  #         * `òperation$sizes`: A vector of integers specifying the sizes of the new groups, being
  #         the first element the remaining size of the original group.
  #  * Drift between groups: It moves a fixed number of agents between the groups
  #         specified and in the order in which they appear. For instance it could
  #         move to elements from group 1 to 2, another 2 (different ones)  from 2 to 3, 
  #         and so on.
  # * `operation$type = "drift"`. 
  # * `operation$groups`: An ordered vector of integers specifying the groups transfering their elements.
  # * `operation$size`: An integer indicating the number of elements transferred between groups.
  test1 = unique(assignments)
  test2 = unique(operation$groups)
  matched = match(test2, test1)
  if(any(is.na(matched)) == T){
    stop(">> groups: ", test2[is.na(matched)], " do not exist in operation ",operation)
  }

  if (operation$type == "merge") {
    # all groups specified will be merged in the first  group in the list, there is no relabelling
    # for the remainder groups
    for (i in 2:length(operation$groups)) {
      assignments[assignments == operation$groups[i]] <- operation$groups[1]
    }
  } else if (operation$type == "split") {
    subset_indices <- which(assignments == operation$group)
    if (length(subset_indices) != sum(operation$sizes)) {
      stop("Split sizes do not match subset size.")
    }
    # A new id is created for each subgroup, the last subgroup keeps the id 
    #browser()
    nmax = max(assignments)
    for(k in 1:(length(operation$sizes)-1)){
      assignments[subset_indices[1:operation$sizes[k]]] <- nmax + k
    }
    #assignments[subset_indices[(operation$sizes[1] + 1):length(subset_indices)]] <- max(assignments) + 2
  } else if (operation$type == "drift") {
    temp <- assignments
    #browser()
    for (i in seq_along(operation$groups)) {
      from <- operation$groups[i]
      to <- operation$groups[(i %% length(operation$groups)) + 1]
      move_indices <- sample(which(assignments == from), operation$size)
      temp[move_indices] <- to
    }
    assignments <- temp
  }
  return(assignments)
}



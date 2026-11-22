# In this problem there is one merge at each step
#problem = "merge"
N <- 200  # Number of agents
M <- 10  # Number of probability matrices
K <- 10   # Number of subgroups

subset_sizes <- c(5, 20, 10, 40, 20, 7, 15, 55, 10, 18)  # Sizes of each subset

operations_list = list()
i = 1

#1 groups 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
operations_list[[i]] <- list(
  list(type="merge", groups=c(3, 5))) # 5 dissapears
i = i + 1

#2 groups 1, 2, 3, 4, 6, 7, 8, 9, 10
operations_list[[i]] <- list(
  list(type="merge", groups=c(8, 1))) # 1 dissapears
i = i + 1

#3 groups 2, 3, 4, 6, 7, 8, 9, 10
operations_list[[i]] <- list(
  list(type="merge", groups=c(8, 6))) # 6 dissapears
i = i + 1

#4 groups 2, 3, 4, 7, 8, 9, 10
operations_list[[i]] <- list(
  list(type="merge", groups=c(4, 10))) # 10 dissapears
i = i + 1

#5 groups 2, 3, 4, 7, 8, 9
operations_list[[i]] <- list(
  list(type="merge", groups=c(2, 9))) # 9 dissapears
i = i + 1

#6 groups 2, 3, 4, 7, 8
operations_list[[i]] <- list(
  list(type="merge", groups=c(4, 7))) # 7 dissapears
i = i + 1

#7 groups 2, 3, 4, 8
operations_list[[i]] <- list(
  list(type="merge", groups=c(3, 4))) # 4 dissapears
i = i + 1

#8 groups 2, 3, 8
operations_list[[i]] <- list(
  list(type="merge", groups=c(2, 8))) # 8 dissapears
i = i + 1

#9 groups 2, 3
operations_list[[i]] <- list(
  list(type="merge", groups=c(2, 3))) # 3 dissapears
i = i + 1
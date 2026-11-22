# In this problem there is one merge at each step
#problem = "merge-split-drift"
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

#2 groups 1, 2, 3, 4, 6, 7, 8, 9, 10, 11
operations_list[[i]] <- list(
  list(type="split", group=4, sizes = c(10,30))) # 11 appears
i = i + 1

#3 groups 1, 2, 3, 4, 6, 7, 8, 9, 10, 11
operations_list[[i]] <- list(
  list(type="drift", groups=c(1, 3, 6, 8, 10), size=2)) # drift
i = i + 1

#4 groups 1, 2, 3, 4, 6, 7, 8, 9, 10, 11
operations_list[[i]] <- list(
  list(type="merge", groups=c(8, 6))) # 6 dissapears
i = i + 1

#5 groups 1, 2, 3, 4, 7, 8, 9, 10, 11
operations_list[[i]] <- list(
  list(type="merge", groups=c(4, 10))) # 10 dissapears
i = i + 1

#6 groups 1, 2, 3, 4, 7, 8, 9, 11
operations_list[[i]] <- list(
  list(type="drift", groups=c(2, 4, 8, 11, 1), size=4)) # drift
i = i + 1

#7 groups 1, 2, 3, 4, 7, 8, 9, 11, 12, 13, 14 
operations_list[[i]] <- list(
  list(type="split", group = 8, sizes = c(12,22,20,8))) # 12, 13, 14, appear 
i = i + 1

#8 groups 1, 2, 3, 4, 7, 8, 9, 11, 12, 13, 14 
operations_list[[i]] <- list(
  list(type="merge", groups=c(1,2))) # 2 dissapear
i = i + 1

#9 groups 1, 3, 4, 8, 9, 11, 12, 13, 14 
operations_list[[i]] <- list(
  list(type="merge", groups=c(3, 4, 7, 8))) # 4, 7, 8 dissapear
i = i + 1

# #8 groups 1, 3, 9, 11, 12, 13, 14
# operations_list[[i]] <- list(
#   list(type="merge", groups=c(3, 9))) # 9 dissapears
# i = i + 1
# 
# #9 groups 1, 3, 11, 12, 13, 14
# operations_list[[i]] <- list(
#   list(type="merge", groups=c(1, 3))) # 3 dissapears
# i = i + 1


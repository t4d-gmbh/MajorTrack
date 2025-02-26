rm(list = ls())
########### START EDITING
problem = "problem_merge"
p_intra = 0.2
p_inter = 0.05
time_total = 100
time_interv = 10
time_steps <- rep(time_interv, time_total/time_interv)  # Number of time steps for each probability matrix
seed = 25022025 # today
########### STOP EDITING

set.seed(seed)

# Set up environment and source functions --------
library(rprojroot)
dir.project = find_rstudio_root_file()
setwd(dir.project)

dir.src = "src/"
dir.data = "data/"

setwd(dir.src)
problem_code = paste0(problem,".R")
library(gplots)
source("generate_modify_partitions.R")
source("generate_probability_matrices.R")
source("simulate_interactions.R")
source(problem_code)
setwd(dir.project)

setwd(dir.data)
label = paste0(problem,'_N-',N,"_K-",K,"_pin-",p_intra,"_pbet-",p_inter)
dir.create(label)
setwd(label)

# Create structures ----------
# --- Create assignments 
# ...... Initialize assignments randomly given the specified structure
assignments_list <- lapply(1:M, function(x) generate_partitions(N, K, subset_sizes))

# ..... Modify assignments, considering the first assignment, modify the following ones according to operations list
for (m in 2:M) {
  operations_tmp = operations_list[[m - 1]]
  for(u in 1:length(operations_tmp)){
    if(u == 1){ # the first operation should consider the assignment in previous step
      step = m - 1
    }else{ # the remainder, the present  step
      step = m
    }
    assignments_list[[m]] <- modify_partitions(assignments_list[[step]], operations_tmp[[u]])
  }
}

# Create probability matrices ----
# --- Determine structure of probability matrices
Kend = max(unlist(lapply(assignments_list,FUN = max))) #   assignments_list[[length(assignments_list)]])
p_intra_vec = rep(p_intra, Kend)
p_inter_mat = matrix(rep(p_inter, Kend^2), nrow=Kend, byrow=TRUE)
diag(p_inter_mat) = 0

p_intra_list <- replicate(M, p_intra_vec, simplify=FALSE)  # Different intra-group probabilities for each matrix
p_inter_list <- replicate(M, p_inter_mat, simplify=FALSE)  

prob_matrices <- generate_probability_matrices(N, M, assignments_list, p_intra_list, p_inter_list)

# Simulate --------
interaction_matrices <- simulate_interactions(N, prob_matrices, time_steps)


# Double check -------

err_mat = list()
err_tot = c()
t  = 0
for(i in 1:(M-1)){
  int_tmp = 0
  for(k in 1:time_steps[i]){
    t  = t + 1
    int_tmp = int_tmp + interaction_matrices[[t]]
    #file_int_mat = paste0("matrix_int-",t,"_",label,".txt")
    #write.table(interaction_matrices[[t]],
    #            file = file_int_mat, sep = "\t")
  }
  int_tmp = int_tmp / time_steps[i]
  err_mat[[i]] = prob_matrices[[i]] - int_tmp
  err_tot[i] = sum(err_mat[[i]]^2)/(N^2)
  file_prob = paste0("heatmap_prob-",i,"_",label,".pdf")
  file_mat = paste0("heatmap_mat-",i,"_",label,".pdf")
  pdf(file_prob)
  heatmap.2(prob_matrices[[i]],
            dendrogram = "none",
            scale = "none",
            trace = "none",
            key.title = "", 
            #keysize = 0.07, 
            density.info = "none",
            cexRow = 0.2, cexCol = 0.2)
  dev.off()
  pdf(file_mat)
  heatmap.2(interaction_matrices[[t]],
            dendrogram = "none",
            scale = "none",
            trace = "none",
            key.title = "", 
            #keysize = 0.07, 
            density.info = "none",
            cexRow = 0.2, cexCol = 0.2)
  dev.off()
  
}

# Final plots and save data ---------
file_err = paste0("Err_",label,".pdf")
pdf(file_err)
plot(seq(1, length(err_tot)),err_tot,xlab = "Prob. matrix id.")
dev.off()

file_int_mat = paste0("matrices_interaction_",label,".RDS")
saveRDS(interaction_matrices, file_int_mat)

file_assign = paste0("assignments_",label,".RDS")
saveRDS(assignments_list, file_assign)

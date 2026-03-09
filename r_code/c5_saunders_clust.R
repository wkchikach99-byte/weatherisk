args = commandArgs(trailingOnly=TRUE)

path_out <- paste0("~/maxstable/complete/",args[1])

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))

a_sim_exp_ns <- readRDS(
  paste0(path_out,"/res/simulation.rds"))

v_matrix <- c_extrcoeff_matrix(a_sim_exp_ns, madogram=T)
saveRDS(v_matrix, paste0(path_out, "/res/distmatrix_saunders.rds"))

hc <- clustering(v_matrix)
# k_clust <- cluster_number_threshold_method(hc, 0.15)
v_matrix_tri <- v_matrix[upper.tri(v_matrix)]
k_clust <- cluster_number_threshold_method(hc,
                          quantile(v_matrix_tri, 0.3))
saveRDS(k_clust, paste0(path_out, "/res/kclust_saunders.rds"))
clusters <- cutree(hc, k_clust)
saveRDS(clusters, paste0(path_out, "/res/clusters_saunders.rds"))
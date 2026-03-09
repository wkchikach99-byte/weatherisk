args = commandArgs(trailingOnly=TRUE)

path_out <- paste0("~/maxstable/complete/",args[1])
df <- as.numeric(args[2])
alpha   <- as.numeric(args[3])/10

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))

local_estimates_matrix_smoothed <- readRDS(
  paste0(path_out,"/res/local_estimates_",args[2],"_",args[3],"_sm.rds"))

ec_distance_matrix_s <- calc_distance_ellipses(
  local_estimates_matrix_smoothed)
 saveRDS(ec_distance_matrix_s, paste0(path_out, "/res/distmatrix_locest_",args[2],"_",args[3],".rds"))
hc <- clustering(ec_distance_matrix_s)
# k_clust <- cluster_number_threshold_method(hc, 20)
# alternatively: fixed number of clusters:

distance_sim <- function(lat1, lon1, lat2, lon2) {
  return(sqrt((lat1-lat2)^2 + (lon1-lon2)^2))
}
ec_distance_matrix_s_tri <- ec_distance_matrix_s[upper.tri(ec_distance_matrix_s)]
k_clust <- cluster_number_threshold_method(hc,
                    quantile(ec_distance_matrix_s_tri, 0.3))

saveRDS(k_clust, paste0(path_out, "/res/kclust_locest_",args[2],"_",args[3],".rds"))
clusters <- cutree(hc, k_clust)
saveRDS(clusters, paste0(path_out, "/res/clusters_locest_",args[2],"_",args[3],".rds"))
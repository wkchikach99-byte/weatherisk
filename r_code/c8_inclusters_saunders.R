library(lhs)

args = commandArgs(trailingOnly=TRUE)

path_out <- paste0("~/maxstable/complete/",args[1])
df <- as.numeric(args[2])
alpha   <- as.numeric(args[3])/10
clusternum <- as.numeric(args[4])

##################################### Calculations #########################################

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))


a_sim_exp_ns <- readRDS(paste0(path_out,"/res/simulation.rds"))
k_clust <- readRDS(paste0(path_out,"/res/kclust_saunders.rds"))
print(paste("Cluster", clusternum, "of", k_clust))

if(clusternum <= k_clust)
{
	if(alpha < 1) upper_bounds <- c(max(a_matrix)+5, 4*max(b_matrix))
	if(alpha ==1) upper_bounds <- c(max(a_matrix)+5, 2*max(b_matrix))
	if(alpha > 1) upper_bounds <- c(max(a_matrix)+5, 1*max(b_matrix))

	cluster_file = "clusters_saunders"
	print(Sys.time())
	cat(paste0("File '",cluster_file,"'. df = ",df,", alpha = ",alpha, "\n"))
	cluster <- readRDS(paste0(path_out,"/res/",cluster_file,".rds"))
	res <- calc_estimates_in_clusters(cluster, df, alpha, upper_bounds,
									  clusternum=clusternum)
	cat("Saving file.\n")
	saveRDS(res, paste0(path_out,"/res/tmp/in",cluster_file,"_",args[2],"_",args[3],"_",
						clusternum,".rds"))
	print(Sys.time())
}
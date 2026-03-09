args = commandArgs(trailingOnly=TRUE)

path_out <- paste0("~/maxstable/complete/",args[1])
df <- as.numeric(args[2])
alpha   <- as.numeric(args[3])/10

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))

path_in = paste0(path_out,"/res/tmp/")
k_clust <- readRDS(paste0(path_out,"/res/kclust_locest_",args[2],"_",args[3],".rds"))

tmp <- t(matrix(-Inf, nrow=5, ncol=k_clust))
for(k in 1:k_clust) {
  tmp2 <- readRDS(paste0(path_in,"inclusters_locest_",args[2],"_",args[3],"_",k,".rds"))
  tmp <- pmax(tmp,tmp2)
}
tmp[which(is.infinite(tmp))] <- 0
saveRDS(tmp, paste0(path_out,"/res/inclusters_locest_",args[2],"_",args[3],".rds"))
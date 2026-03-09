args = commandArgs(trailingOnly=TRUE)

path_out <- paste0("~/maxstable/complete/",args[1])
df <- as.numeric(args[2])
alpha   <- as.numeric(args[3])/10

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))

if(length(args) > 3)
{
  smoothing_dist = as.numeric(args[4])
}

path_in = paste0(path_out,"/res/tmp/")
filelist <- list.files(path_in, pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
filelist <- filelist[which(startsWith(filelist,
              paste0("total_estimates_",args[2],"_",args[3],"_")))]
print("Filelist:")
print(filelist)
if(length(filelist)==0) {stop(paste0("Error: No files found under ",
       filelist, "total_estimates_",args[2],"_",args[3],"_"))}

for(i in 1:length(filelist)) {
  assign(paste0("t",i), readRDS(paste0(path_in,filelist[i])))
}

local_estimates_matrix <- matrix(-Inf, nrow=nrow(t1),ncol=ncol(t1))
for(i in 1:nrow(t1)) for(j in 1:ncol(t1)) {
  for(k in 1:(length(filelist))) {
    local_estimates_matrix[i,j] <- max(local_estimates_matrix[i,j],
                                       eval(parse(text=paste0("t",k)))[i,j],na.rm=T)
  }
}
print(paste0("Range of resulting file:", min(local_estimates_matrix), ", ",
             max(local_estimates_matrix)))
saveRDS(local_estimates_matrix, paste0(path_out,"/res/local_estimates_",args[2],"_",args[3],".rds"))
if(is.infinite(min(local_estimates_matrix))) {stop("Missing values detected!")}

print(paste0("Smoothing with distance:", smoothing_dist))
local_estimates_matrix_smoothed <- smooth_local_estimates(
  local_estimates_matrix, smoothing_dist=smoothing_dist)
saveRDS(local_estimates_matrix_smoothed, paste0(path_out,
                                                "/res/local_estimates_",args[2],"_",args[3],"_sm.rds"))
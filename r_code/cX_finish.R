args = commandArgs(trailingOnly=TRUE)

path_out <- paste0("~/maxstable/complete/",args[1])

source(paste0("~/maxstable/complete/parameters.R"))
source(paste0("~/maxstable/complete/functions.R"))

path_in = paste0(path_out,"/res/")
filelist <- list.files(path_in, pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
filelist <- filelist[which(startsWith(filelist,
              paste0("inclusters_")))]
print("Filelist:")
print(filelist)
if(length(filelist)==0) {print(paste0("Path: ", path_in));
  stop(paste0("Error: No files found under ",filelist, "inclusters_"))}

alpha_list <- df_list <- numeric()
for(i in 1:length(filelist)) {
  df_list <-    c(df_list,    strsplit(filelist[i],"_")[[1]][3])
  alpha_list <- c(alpha_list, substr(strsplit(filelist[i],"_")[[1]][4],1,2))
}
df_list <- unique(df_list)
alpha_list <- unique(alpha_list)

results <- matrix(nrow=length(df_list)*length(alpha_list), ncol = 5)
colnames(results) <- c("df","alpha","saunders","new","diff")
cnt=1; for(df in df_list) for(alpha in alpha_list) {
  results[cnt,1] <- as.numeric(df)
  results[cnt,2] <- as.numeric(alpha)/10
  
  tmp <- readRDS(paste0(path_in,"inclusters_saunders_",df,"_",alpha,".rds"))
  wh <- which(tmp[,5]!=0)
  results[cnt,3] <- sum(tmp[wh,4]*tmp[wh,5])/sum(tmp[wh,4])
  
  tmp <- readRDS(paste0(path_in,"inclusters_locest_",df,"_",alpha,".rds"))
  wh <- which(tmp[,5]!=0)
  results[cnt,4] <- sum(tmp[wh,4]*tmp[wh,5])/sum(tmp[wh,4])
  
  results[cnt,5] <- -results[cnt,3]+results[cnt,4]
  cnt = cnt+1
}

filelist <- list.files(path_in, pattern="*.rds", all.files=FALSE,
                       full.names=FALSE)
for(f in filelist) {
  assign(strsplit(f,"\\.")[[1]][1], readRDS(paste0(path_in,f)))
}
rm(f)


saveRDS(results, paste0(path_out, "/results.rds"))
write.table(results, file=paste0(path_out,"/results.txt"), row.names=FALSE, sep="\t")
save.image(file=paste0(path_out,"/",args[1],"_tmp.rdata"))


##################################

calculate_eici <- function(df_num, alpha_num) {
  eici <- matrix(-Inf, nrow=max(clusters1)*max(clusters2), ncol=5)
  cnt=0
  for(i in 1:max(clusters1)) {
    for(j in 1:max(clusters2)) {
      which_cl <- which(clusters1==i & clusters2==j)
      cnt=cnt+1
      print(c(i,j, length(which_cl)))
      if(length(which_cl)>=5) {
        eici[cnt,] <- c(i,j,length(which_cl),
                        llh_in_cluster(
                          t(sapply(which_cl, function(k) {
                            simulation[number_grid(k)[1],number_grid(k)[2],]})),
                          df_num, alpha_num, X[which_cl], Y[which_cl] , 
                          incl1[i,1:3],
                          max_dist=4*((max(X)-min(X))/(resolution-1))),
                        llh_in_cluster(
                          t(sapply(which_cl, function(k) {
                            simulation[number_grid(k)[1],number_grid(k)[2],]})),
                          df_num, alpha_num, X[which_cl], Y[which_cl] , 
                          incl2[j,1:3],
                          max_dist=4*((max(X)-min(X))/(resolution-1))))
      }
    }
  }
  return(eici)
}
for(df in df_list) for(alpha in alpha_list) {
  clusters1 <- clusters_saunders
  clusters2 <- get(paste0("clusters_locest_",df,"_",alpha))
  incl1 <- get(paste0("inclusters_saunders_",df,"_",alpha))
  incl2 <- get(paste0("inclusters_locest_",df,"_",alpha))
  eici <- calculate_eici(as.numeric(df),as.numeric(alpha)/10 )
  saveRDS(eici, paste0(path_in, "eici_",df,"_",alpha,".rds"))
}
filelist <- list.files(path_in, pattern="*.rds", all.files=FALSE,
                       full.names=FALSE)
for(f in filelist) {
  assign(strsplit(f,"\\.")[[1]][1], readRDS(paste0(path_in,f)))
}
rm(f)
save.image(file=paste0(path_out,"/",args[1],"_final.rdata"))
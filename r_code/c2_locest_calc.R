library(lhs)

args = commandArgs(trailingOnly=TRUE)

start <- end <- as.integer(args[1])
path_out <- paste0("~/maxstable/complete/",args[2])
df <- as.numeric(args[3])
alpha   <- as.numeric(args[4])/10

##################################### Calculations #########################################
source("~/maxstable/complete/parameters.R")
source("~/maxstable/complete/functions.R")

# override value in parameters.R
if(length(args) > 4)
{
  locest_abst = as.numeric(args[5])
}

if(start <= resolution) {
  if(alpha < 1) upper_bounds <- c(max(a_matrix)+5, 6*max(b_matrix))
  if(alpha ==1) upper_bounds <- c(max(a_matrix)+5, 2*max(b_matrix))
  if(alpha > 1) upper_bounds <- c(max(a_matrix)+5, 1*max(b_matrix))
  
  dat <- readRDS(paste0(path_out,"/res/simulation.rds"))
  
  local_estimates_matrix <- matrix(NA, nrow=length(X),ncol=3)
  print(Sys.time())
  cat(paste0("Calculating Local Estimates (From: ",start,". To: ",end,")\n"))
  cat(paste0("locest_abst =  ",locest_abst,"\n"))
  
  for(i in (start-1)*resolution+(1:resolution)) {
    local_estimates_matrix[i,] <- tryCatch(pairwise_density_optim_local(dat,
                                     df,alpha, X[i],Y[i] ,abstand=locest_abst, print=F,
                                     upper_bounds = upper_bounds,ensemble = locest_ensemble),
                                           error=function(cond) {message(cond);return(numeric(3)*NA)})
    toprint <- c(i, Y[i], X[i],local_estimates_matrix[i,])
    cat(format(round(toprint,4),nsmall=4)); cat("\n")
  }
  cat("Saving file.\n")
  saveRDS(local_estimates_matrix, paste0(path_out,
                                         "/res/tmp/total_estimates_",args[3],"_",args[4],"_",start,".rds"))
  print(Sys.time())
}
library(lhs)
library(lattice)
library(RColorBrewer)
library(pals)

path = "C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/3. Max-stable/R/results_final/"

oldstuff <- function() {
name = "new_stripes_e4_m2"

path = paste0("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/",
              "3. Max-stable/R/results_varying_distance_smoothing/", name, "/")
load(paste0(path, name, ".RData"))
all_incl <- numeric(0)
for(df in df_list) for(alpha in alpha_list) {
  incl <- get(paste0("inclusters_locest_",df,"_",alpha))
  all_incl <- c(all_incl, incl[which(incl[,5]!=0),5])
  incl <- get(paste0("inclusters_saunders_",df,"_",alpha))
  all_incl <- c(all_incl, incl[which(incl[,5]!=0),5])
}
print(range(all_incl))
hist(all_incl)

von = -4.25
bis = -4

## plot one simulation
# load("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/3. Max-stable/R/results/stripes_cont_5_2/stripes_cont_5_2.rdata")
# plot_map(exp(-simulation[,,1]^(-1)), von=0, bis=1, cr=rev(brewer.pal(11,"RdBu")))



############### example ec's around a fixed point #######################

# plot real values
png(paste0(path,"nonstat34.png"), width=1300, height=1000, unit="px", res=300)
plot_mado_nonstat(-3,  2 ,df_true,alpha_true, a_matrix,b_matrix,g_matrix)
plot_mado_nonstat(-3, -2 ,df_true,alpha_true, a_matrix,b_matrix,g_matrix)
plot_mado_nonstat( 3,  2 ,df_true,alpha_true, a_matrix,b_matrix,g_matrix)
plot_mado_nonstat( 3, -2 ,df_true,alpha_true, a_matrix,b_matrix,g_matrix)
dev.off()

# plot local estimates
#plot_map(local_estimates_5_10_sm[,3])
#extrcoeff_matrix <- c_extrcoeff_matrix(simulation)
#plot_local_estimates(-2, -2, extrcoeff_matrix, local_estimates_5_10_sm)

plot_inclusters <- function(clusters, inclusters, parnum, von, bis)
{
  plot_map(0*X+sapply(clusters, function(i) {inclusters[i,parnum]}),
           von=von, bis=bis,
           cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
           contours=clusters, truncate=T,legend=F)
}

for(df in df_list) for(alpha in alpha_list) {
  png(paste0(path,"estimates_locest_",df,"_",alpha,".png"),
      width=1000, height=1000, unit="px", res=300)
  matr <- get(paste0("clusters_locest_",df,"_",alpha))
  incl <- get(paste0("inclusters_locest_",df,"_",alpha))
  plot_map(0*X+sapply(matr, function(i) {incl[i,3]}),
           von=von, bis=bis,
           cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
           contours=matr, truncate=T,legend=F)
  dev.off()
}

matr <- clusters_saunders
incl <- get(paste0("inclusters_saunders_",df_true,"_",alpha_true*10))

png(paste0(path,"estimates_saunders_",df_true,"_",alpha_true*10,".png"),
    width=1200, height=1000, unit="px", res=300)
plot_map(0*X+sapply(matr,function(i) {incl[i,5]}),
         von=von, bis=bis,
         cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
         contours=matr, truncate=T)
dev.off()

############################ LLH on intersections of clusters #####################################

clusters1 <- clusters_locest_5_10
clusters2 <- clusters_saunders
incl1 <- inclusters_locest_5_10
incl2 <- inclusters_saunders_5_10

clusters_combined <- 10*clusters1 + clusters2

# estimates_in_clusters_intersect
eici <- matrix(-Inf, nrow=max(clusters1)^2, ncol=5)
cnt=0
for(i in 1:max(clusters1)) {
  for(j in 1:max(clusters2)) {
    print(Sys.time())
    print(c(i,j))
    which_cl <- which(clusters1==i & clusters2==j)
    cnt=cnt+1
    if(length(which_cl)>0) {
      eici[cnt,] <- c(i,j,length(which_cl),
        llh_in_cluster(
          t(sapply(which_cl, function(k) {
            simulation[number_grid(k)[1],number_grid(k)[2],]})),
          df_true, alpha_true, X[which_cl], Y[which_cl] , 
          incl1[i,1:3],
          max_dist=locest_abst*((max(X)-min(X))/(resolution-1))),
        llh_in_cluster(
          t(sapply(which_cl, function(k) {
            simulation[number_grid(k)[1],number_grid(k)[2],]})),
          df_true, alpha_true, X[which_cl], Y[which_cl] , 
          incl2[j,1:3],
          max_dist=locest_abst*((max(X)-min(X))/(resolution-1))))
    }
  }
}
}

############### threshold to euclidean distance #################
distance_geo <- function(lat1, lon1, lat2, lon2) {
  return(acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos((lon2-lon1)*pi/180))*6371)
}
# or for simulations:
distance_sim <- function(lat1, lon1, lat2, lon2) {
  return(sqrt((lat1-lat2)^2 + (lon1-lon2)^2))
}

plot_threshold_distances <- function(distmatrix, x_max, y_max, dist_fkt, yt) {
  cr=colorRampPalette(brewer.pal(11,"RdBu")[6:11])(21)
  heatmap_xax <- seq(0.2, x_max, length.out = 35)
  heatmap_yax <- seq(0, y_max, length.out=100)
  heatmap <- matrix(0, nrow = length(heatmap_xax), ncol = length(heatmap_yax))
  
  for(i in 2:length(X)) for (j in 1:(i-1))
  {
    dist = dist_fkt(Y[i], X[i], Y[j],X[j])
    x_pos = length(which(dist - heatmap_xax >= 0))  # at which position in heatmap_xax is the value dist?
    res = distmatrix[i,j]
    y_pos = length(which(res - heatmap_yax >= 0))
    if(res <= y_max + heatmap_yax[2] - heatmap_yax[1])
    {
      if(dist <= x_max + heatmap_xax[2] - heatmap_xax[1])
      {
        heatmap[x_pos, y_pos] = heatmap[x_pos, y_pos] + 1
      }
    }
  }
  
  # normalize number of values
  # for(i in 1:dim(heatmap)[1])
  # {
  #   #heatmap[i,] = heatmap[i,]/sum(heatmap[i,])
  #   #heatmap[i,] <- pmin(0.15, heatmap[i,])
  # }
  
  d = data.frame(x=rep(heatmap_xax, ncol(heatmap)), 
                 y=rep(heatmap_yax, each=nrow(heatmap)), 
                 z=heatmap)
  levelplot(heatmap~x*y, data = d, col.regions = rev(cr),
            xlab="Euclidean Distance", ylab= yt)
}

# simulations:
png(paste0(path,"/thr_euclid_locest.png"),
    width=1200, height=1000, unit="px", res=300)
plot_threshold_distances(distmatrix_saunders, 4, 0.2, distance_sim, "Dissimilarity measure")
dev.off()
png(paste0(path, "/thr_euclid_saunders.png"),
    width=1200, height=1000, unit="px", res=300)
plot_threshold_distances(distmatrix_locest_5_10, 4, 80, distance_sim, "Dissimilarity measure")
dev.off()
# Fig 3


## quantiles
ec_distance_matrix_s_tri <- distmatrix_saunders[upper.tri(distmatrix_saunders)]
quantile(ec_distance_matrix_s_tri, 0.3)
ec_distance_matrix_s_tri <- distmatrix_locest_5_10[upper.tri(distmatrix_locest_5_10)]
quantile(ec_distance_matrix_s_tri, 0.3)

##################################### inclusters ########################################

### stripes:
# nm = "stripes_"
# index = 2
# von = 0; bis = 5
# m = b_matrix
### bigsmall:
# nm = "bigsmall_"
# index = 1
# von = 1; bis = 5
# m = a_matrix
### rotate:
nm = "rotate_"
index = 3
von = -pi/2; bis = pi/2
m = g_matrix
###

alpha_true = alpha; df_true = df
png(paste0(path,nm,"estimates_locest_",df_true,"_",alpha_true,".png"),
    width=1000, height=1000, unit="px", res=300)
matr <- get(paste0("clusters_locest_",df_true,"_",alpha_true))
incl <- get(paste0("inclusters_locest_",df_true,"_",alpha_true))
plot_map(0*X+sapply(matr, function(i) {incl[i,index]}),
         von=von, bis=bis,
         cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
         contours=matr, truncate=T,legend=F)
dev.off()

matr <- clusters_saunders
incl <- get(paste0("inclusters_saunders_",df_true,"_",alpha_true))

png(paste0(path,nm,"estimates_saunders_",df_true,"_",alpha_true,".png"),
    width=1000, height=1000, unit="px", res=300)
plot_map(0*X+sapply(matr,function(i) {incl[i,index]}),
         von=von, bis=bis,
         cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
         contours=matr, truncate=T,legend=F)
dev.off()

# real values
png(paste0(path,nm,"realvalues.png"),
width=1200, height=1000, unit="px", res=300)
plot_map(m, von=von, bis=bis,cr=colorRampPalette(brewer.pal(11,"RdBu")[6:11])(101))
dev.off()

# Fig 4

################################# betterworse ###########################################

plot_parameters <- function() {
  a_min <- min(a_matrix); a_max <- max(a_matrix)
  b_min <- min(b_matrix); b_max <- max(b_matrix)
  g_min <- min(g_matrix); g_max <- max(g_matrix)
  if(a_min==a_max) {a_min <- a_min-0.5; a_max <- a_max+0.5}
  if(b_min==b_max) {b_min <- b_min-0.5; b_max <- b_max+0.5}
  if(g_min==g_max) {g_min <- g_min-0.5; g_max <- g_max+0.5}
  plot_map(a_matrix, von=a_min, bis=a_max, cr=colorRampPalette(brewer.pal(9,"Blues"))(101))
  #plot_map(b_matrix, von=b_min, bis=b_max, cr=colorRampPalette(brewer.pal(9,"Blues"))(101))
  plot_map(g_matrix, cr=colorRampPalette(brewer.pal(9,"Blues"))(101))
  # parameters saunders
  plot_map(0*X+sapply(clusters1,function(i) {incl1[i,1]}),
           von=a_min, bis=a_max,
           cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
           contours=clusters1, truncate=T)
  #plot_map(0*X+sapply(clusters1,function(i) {incl1[i,2]}),
  #         von=b_min, bis=b_max,
  #         cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
  #         contours=clusters1, truncate=T)
  plot_map(0*X+sapply(clusters1,function(i) {incl1[i,3]}),
           von=-pi/2, bis=pi/2,
           cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
           contours=clusters1, truncate=T)
  # parameters locest
  plot_map(0*X+sapply(clusters2,function(i) {incl2[i,1]}),
           von=a_min, bis=a_max,
           cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
           contours=clusters2, truncate=T)
  #plot_map(0*X+sapply(clusters2,function(i) {incl2[i,2]}),
  #         von=b_min, bis=b_max,
  #         cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
  #         contours=clusters2, truncate=T)
  plot_map(0*X+sapply(clusters2,function(i) {incl2[i,3]}),
           von=-pi/2, bis=pi/2,
           cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
           contours=clusters2, truncate=T)
}

calculate_eici <- function() {
  eici <- matrix(-Inf, nrow=max(clusters1)^2, ncol=5)
  cnt=0
  for(i in 1:max(clusters1)) {
    for(j in 1:max(clusters2)) {
      which_cl <- which(clusters1==i & clusters2==j)
      cnt=cnt+1
      #print(c(i,j, length(which_cl)))
      if(length(which_cl)>0) {
        eici[cnt,] <- c(i,j,length(which_cl),
                        llh_in_cluster(
                          t(sapply(which_cl, function(k) {
                            simulation[number_grid(k)[1],number_grid(k)[2],]})),
                          df_num, alpha_num, X[which_cl], Y[which_cl] , 
                          incl1[i,1:3],
                          max_dist=locest_abst*((max(X)-min(X))/(resolution-1))),
                        llh_in_cluster(
                          t(sapply(which_cl, function(k) {
                            simulation[number_grid(k)[1],number_grid(k)[2],]})),
                          df_num, alpha_num, X[which_cl], Y[which_cl] , 
                          incl2[j,1:3],
                          max_dist=locest_abst*((max(X)-min(X))/(resolution-1))))
      }
    }
  }
  return(eici)
}

# results are searched for in 3. Max-Stable/R/[resultspath]/[name_] if num_iter = 1
# if num_iter > 1, the paths are 3. Max-Stable/R/[resultspath]/[name_][i] instead, for 1 <= i <= num_iter

resultspath = "results_final"
name_ = "stripes_cont_av"
num_iter = 25


#name_ = "cont_stripes_all_"
#name_ = "rotate_inv"
#name_ = "bigsmall_"


#for(df_txt in c("3","5","7")) for(alpha_txt in c("07","10","13"))
#{
  df_txt = "5"; alpha_txt = "10"
  df_num <- as.numeric(df_txt)
  alpha_num <- as.numeric(alpha_txt)/10
  #print(Sys.time())
  print(paste0("df = ",df_num,", alpha = ",alpha_num))
  comp <- matrix(nrow=num_iter, ncol=2)
  bw = 0*X
  for(n_ in 1:num_iter) {
    print(paste0("   ",Sys.time()))
    if(num_iter > 1)
    {
      print(paste0("   ",n_))
      name = paste0(name_, n_)
    }
    else
    {
      name = name_
    }
    path = paste0("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/",
                  "3. Max-stable/R/", resultspath, "/", name, "/")
    load(paste0(path, name, "_final.RData"))
    source(paste0("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/",
                  "3. Max-stable/R/codefiles/functions.R"))
    
    eicimap <- function(clusters_saunders, clusters_locest, eici)
    {
      clusters_combined <- 1000*clusters_saunders + clusters_locest
      return(0*X+sapply(clusters_combined, function(i) {
        k=which(eici[,2]==i%%1000 & eici[,1]==floor(i/1000))[1];
        (eici[k,5]>eici[k,4])
      }))
    }
    
    eici <- get(paste0("eici_",df_txt,"_",alpha_txt))
    clusters1 <- clusters_saunders
    clusters2 <- get(paste0("clusters_locest_",df_txt,"_",alpha_txt))
    clusters_combined <- 1000*clusters1 + clusters2
    incl1 <- get(paste0("inclusters_saunders_",df_txt,"_",alpha_txt))
    incl2 <- get(paste0("inclusters_locest_",df_txt,"_",alpha_txt))
    
    png(paste0(path,"clusters1_",df_txt,"_",alpha_txt,".png"),
            width=1000, height=1000, unit="px", res=300)
    plot_map(clusters1, von=1, bis=max(clusters1), contours=clusters1,legend=F)
    dev.off()
    
    png(paste0(path,"clusters2_",df_txt,"_",alpha_txt,".png"),
            width=1000, height=1000, unit="px", res=300)
    plot_map(clusters2, von=1, bis=max(clusters2), contours=clusters2,legend=F)
    dev.off()
    
    # plot parameters
    #plot_parameters()
    
    #eici <- readRDS(paste0(path,paste0(path,"eici_",df_txt,"_",alpha_txt,".rds")))
     png(paste0(path,"bw_",df_txt,"_",alpha_txt,".png"),
         width=1000, height=1000, unit="px", res=300)
    eicimap <- eicimap(clusters1, clusters2, eici)
    plot_map(eicimap, von=-0.5, bis=+1.5,
             cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
             contours=clusters_combined, truncate=F,legend=F)
    dev.off() # Fig5
    comp[n_,1] <- sum(eici[which(is.finite(eici[,4])),4])
    comp[n_,2] <- sum(eici[which(is.finite(eici[,5])),5])
    bw <- bw + eicimap/num_iter
  }
  if (num_iter > 1)
  {
    print(colSums(comp))
    png(paste0("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/",
             "3. Max-stable/R/", resultspath, "/",name,"/", "betterworse",df_txt,"_",alpha_txt,".png"),
             width=1200, height=1000, unit="px", res=300)
    plot_map(bw*100, von=0, bis=100, cr=colorRampPalette(brewer.pal(9,"Blues"))(101))
    dev.off() # Fig6
  }
#}

  
  ############################# check influence of locest and smoothing distance ###########################
  alpha_true = 10
  for (e in c(2,3,4,5,6))
    {
    name_ = paste0("e_param_",e)
    path = paste0("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/",
                  "3. Max-stable/R/results_final/", name_, "/")
    if(file.exists(paste0(path, name_, "_final.RData"))) {
      load(paste0(path, name_, "_final.RData"))
      alpha_true = 10
      matr <- get(paste0("clusters_locest_",df_true,"_",alpha_true))
      incl <- get(paste0("inclusters_locest_",df_true,"_",alpha_true))
      png(paste0("C:/Users/Justus/Desktop/Doktorarbeit/Doktorarbeit/",
                 "3. Max-stable/R/results_final/", "eparam_",e,".png"),
          width=1000, height=1000, unit="px", res=300)
      plot_map(0*X+sapply(matr, function(i) {incl[i,2]}),
             von=0, bis=5,
             cr=colorRampPalette(brewer.pal(9,"Blues"))(101),rev=T,
             contours=matr, truncate=T,legend=F)
      dev.off()
      #plot_cluster_map(clusters_locest_5_10, von=1, bis=5, contours=clusters_locest_5_10,legend=F, main=name_)
    }
  }

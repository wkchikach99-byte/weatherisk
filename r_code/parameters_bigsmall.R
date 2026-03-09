resolution = 51 # also hard-coded in c2_locest_calc.sh !
df_true = 5
alpha_true = 1
n_sim = 250


x_ax <-  seq(-5,5, length.out=resolution)   # must be from low to high
y_ax <- -seq(-5,5, length.out=resolution)  # must be from high to low

X <- rep(1,length(y_ax))%*%t(x_ax)
Y <- y_ax%*%t(rep(1,length(x_ax)))

nrow_grid <- length(y_ax)
ncol_grid <- length(x_ax)
n_grid <- nrow_grid*ncol_grid
exp.grid <- expand.grid(y=y_ax,x=x_ax)

### matrix parameters ###


a_matrix <- (7.5 - sqrt((X^2)+Y^2))/2 + 1
b_matrix <- 0*X
g_matrix <- 0*X

locest_ensemble = 5 # number of iterations for optim
locest_abst = 4
smoothing_dist = 2


### testing of parameters:
# a_matrix <- (7.5 - sqrt((X^2)+Y^2))/2+1
# b_matrix <- 0*X
# g_matrix <- 0*X
# 
# plot_map(a_matrix)
# 
# ec_distance_matrix_s <- calc_distance_ellipses(cbind(as.vector(a_matrix),
#                                                      as.vector(b_matrix),
#                                                      as.vector(g_matrix)), res=21)
# hc <- clustering(ec_distance_matrix_s)
# clusters <- cutree(hc, k_clust)
# plot_cluster_map(clusters)

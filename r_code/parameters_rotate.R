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


a_matrix <- 0*X + 1
b_matrix <- 0*X + 3 # from 0 to b_end
g_matrix <- -(X/5 * (pi/2) )

locest_ensemble = 5 # number of iterations for optim
locest_abst = 4
smoothing_dist = 2


## test
# 
# resolution = 11
# df_true = 5
# alpha_true = 1
# n_sim = 10
# 
# 
# x_ax <-  seq(-5,5, length.out=resolution)   # must be from low to high
# y_ax <- -seq(-5,5, length.out=resolution)  # must be from high to low
# 
# X <- rep(1,length(y_ax))%*%t(x_ax)
# Y <- y_ax%*%t(rep(1,length(x_ax)))
# 
# nrow_grid <- length(y_ax)
# ncol_grid <- length(x_ax)
# n_grid <- nrow_grid*ncol_grid
# exp.grid <- expand.grid(y=y_ax,x=x_ax)
# 
# ### matrix parameters ###
# 
# b_end = 5
# 
# a_matrix <- 0*X + 2
# b_matrix <- (X+5)/10*b_end # from 0 to b_end
# g_matrix <- 0*X
# 
# k_clust = 5
# 
# locest_ensemble = 5 # number of iterations for optim
# locest_abst = 4
# smoothing_dist = 2
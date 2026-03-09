#!/usr/bin/env Rscript
###############################################################################
# Generate reference data from the R code for Python verification tests.
# 
# This script runs through every computational step of the R pipeline
# on a small 10x10 grid (stripes preset), saves intermediate results
# as CSV files in tests/reference_data/, and writes a manifest listing
# all outputs with their descriptions.
#
# Usage:  cd climate-risk && Rscript tests/generate_r_reference.R
###############################################################################

cat("=== R Reference Data Generator ===\n")
cat("Working directory:", getwd(), "\n\n")

# ── 0.  Source R code and parameters ────────────────────────────────────────
library(lhs)                      # needed by pairwise_density_optim_local
source("r_code/functions.R")      # all math functions
source("r_code/parameters.R")     # stripes preset (resolution=10, etc.)

outdir <- "tests/reference_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Write manifest header
manifest <- file(file.path(outdir, "MANIFEST.txt"), open = "w")
writeLines("# R Reference Data Manifest", manifest)
writeLines(paste0("# Generated: ", Sys.time()), manifest)
writeLines(paste0("# R version: ", R.version.string), manifest)
writeLines("", manifest)

save_csv <- function(obj, name, description) {
  path <- file.path(outdir, paste0(name, ".csv"))
  write.csv(obj, path, row.names = FALSE)
  writeLines(paste0(name, ".csv  —  ", description), manifest)
  cat("  Saved:", name, "\n")
}

# ── 1.  Grid & Parameters ──────────────────────────────────────────────────
cat("\n[Step 1] Grid & Parameters\n")

# Save grid coordinates (X,Y matrices, column-major flat order)
grid_df <- data.frame(
  flat_index = 1:n_grid,      # R 1-based column-major index
  X = as.vector(X),           # column-major flatten
  Y = as.vector(Y)
)
save_csv(grid_df, "grid_coordinates",
         "Grid coordinates, column-major (R order). Cols: flat_index, X, Y")

# Save parameter matrices (also column-major flat)
params_df <- data.frame(
  flat_index = 1:n_grid,
  a = as.vector(a_matrix),
  b = as.vector(b_matrix),
  g = as.vector(g_matrix)
)
save_csv(params_df, "parameter_matrices",
         "True a, b, gamma parameter matrices (stripes preset), column-major")

# Save scalar parameters
scalars_df <- data.frame(
  parameter = c("resolution", "df_true", "alpha_true", "n_sim",
                 "nrow_grid", "ncol_grid", "n_grid",
                 "locest_ensemble", "locest_abst", "smoothing_dist", "b_end"),
  value = c(resolution, df_true, alpha_true, n_sim,
            nrow_grid, ncol_grid, n_grid,
            locest_ensemble, locest_abst, smoothing_dist, b_end)
)
save_csv(scalars_df, "scalar_parameters",
         "Scalar parameters for stripes preset")

# Save axis vectors  
save_csv(data.frame(x_ax = x_ax), "x_axis", "x-axis values")
save_csv(data.frame(y_ax = y_ax), "y_axis", "y-axis values (high to low)")

# ── 2.  Covariance Function Tests ──────────────────────────────────────────
cat("\n[Step 2] Covariance Function\n")

# Test cov_fkt_2d at specific points
cov_test_cases <- expand.grid(
  x = c(0, 1, -1, 2.5, 0.5),
  y = c(0, 1, -1, 2.5, 0.5),
  a = c(1, 2),
  b = c(0, 1, 3),
  g = c(0, pi/4, -pi/3)
)
cov_test_cases$alpha <- 1.0
cov_test_cases$cov_value <- mapply(
  function(x, y, a, b, g, alpha) cov_fkt_2d(x, y, alpha, a, b, g),
  cov_test_cases$x, cov_test_cases$y,
  cov_test_cases$a, cov_test_cases$b,
  cov_test_cases$g, cov_test_cases$alpha
)
save_csv(cov_test_cases, "cov_fkt_2d_test_cases",
         "cov_fkt_2d evaluated at various (x,y,a,b,g,alpha) combinations")

# Test with alpha != 1
cov_alpha_cases <- data.frame(
  x = c(1, 2, 0.5, 1, 1),
  y = c(1, 0, 1.5, -1, 0.3),
  alpha = c(0.5, 0.5, 1.5, 2, 0.8),
  a = c(1, 2, 1, 3, 1.5),
  b = c(1, 0, 2, 1, 0.5),
  g = c(0, pi/6, -pi/4, pi/3, 0)
)
cov_alpha_cases$cov_value <- mapply(
  function(x, y, alpha, a, b, g) cov_fkt_2d(x, y, alpha, a, b, g),
  cov_alpha_cases$x, cov_alpha_cases$y,
  cov_alpha_cases$alpha, cov_alpha_cases$a,
  cov_alpha_cases$b, cov_alpha_cases$g
)
save_csv(cov_alpha_cases, "cov_fkt_2d_alpha_cases",
         "cov_fkt_2d with various alpha values")

# ── 3.  Non-stationary Covariance ──────────────────────────────────────────
cat("\n[Step 3] Non-stationary Covariance\n")

nonstat_cases <- data.frame(
  x = c(1, 0, 2, -1, 0.5),
  y = c(0, 1, -1, 2, 0.5),
  alpha = c(1, 1, 1, 1, 1),
  a1 = c(2, 1, 3, 1.5, 2),
  b1 = c(1, 2, 0, 1, 3),
  g1 = c(0, pi/4, -pi/6, pi/3, 0),
  a2 = c(1, 2, 2, 2, 1.5),
  b2 = c(3, 1, 1, 0.5, 2),
  g2 = c(pi/6, 0, pi/4, -pi/4, pi/3)
)
nonstat_cases$cov_value <- mapply(
  function(x, y, alpha, a1, b1, g1, a2, b2, g2)
    cov_fkt_2d_nonstat2(x, y, alpha, a1, b1, g1, a2, b2, g2),
  nonstat_cases$x, nonstat_cases$y, nonstat_cases$alpha,
  nonstat_cases$a1, nonstat_cases$b1, nonstat_cases$g1,
  nonstat_cases$a2, nonstat_cases$b2, nonstat_cases$g2
)
save_csv(nonstat_cases, "cov_nonstat_test_cases",
         "cov_fkt_2d_nonstat2 evaluated at various parameter combinations")

# ── 4.  Extremal Coefficient Conversion ────────────────────────────────────
cat("\n[Step 4] Extremal Coefficient Conversion\n")

ec_cases <- expand.grid(
  df = c(1, 3, 5, 10),
  cov = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0)
)
ec_cases$ec <- mapply(cov_to_ec, ec_cases$df, ec_cases$cov)
save_csv(ec_cases, "cov_to_ec_test_cases",
         "cov_to_ec conversion for various (df, cov) combinations")

# ec_to_cov (inverse)
ec_inv_cases <- expand.grid(
  df = c(1, 5, 10),
  ec = c(1.0, 1.1, 1.3, 1.5, 1.7, 1.9)
)
ec_inv_cases$cov <- mapply(ec_to_cov, ec_inv_cases$df, ec_inv_cases$ec)
save_csv(ec_inv_cases, "ec_to_cov_test_cases",
         "ec_to_cov inverse conversion for various (df, ec) combinations")

# ── 5.  Pairwise Density Summand ──────────────────────────────────────────
cat("\n[Step 5] Pairwise Density Summand\n")

# Test density at specific z1, z2, x, y, params
density_cases <- data.frame(
  z1 = c(1.5, 2.0, 5.0, 0.5, 3.0, 10.0, 1.0, 2.5),
  z2 = c(2.0, 1.5, 3.0, 1.0, 4.0, 8.0, 1.5, 3.0),
  x  = c(1.0, 0.5, 2.0, 1.0, -1.0, 0.0, 1.5, 2.0),
  y  = c(0.0, 1.0, 1.0, -1.0, 0.5, 1.0, 0.0, -1.0),
  df = c(5, 5, 5, 5, 5, 3, 10, 5),
  alpha = c(1, 1, 1, 1, 1, 1, 1, 1),
  a  = c(2.0, 1.0, 3.0, 1.5, 2.0, 1.0, 2.0, 2.5),
  b  = c(1.0, 2.0, 0.0, 1.0, 3.0, 1.0, 0.5, 1.5),
  g  = c(0.0, pi/4, 0.0, -pi/6, pi/3, 0, 0, pi/4)
)
density_cases$density_value <- mapply(
  function(z1, z2, x, y, df, alpha, a, b, g)
    pairwise_density_summand(z1, z2, x, y, df, alpha, a, b, g),
  density_cases$z1, density_cases$z2,
  density_cases$x, density_cases$y,
  density_cases$df, density_cases$alpha,
  density_cases$a, density_cases$b, density_cases$g
)
save_csv(density_cases, "pairwise_density_test_cases",
         "pairwise_density_summand for various (z1,z2,x,y,df,alpha,a,b,g)")

# Also test dtdiff
dtdiff_cases <- data.frame(
  x = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
  df = c(6, 6, 6, 6, 6, 6, 6, 6)
)
dtdiff_cases$dtdiff_value <- mapply(dtdiff, dtdiff_cases$x, dtdiff_cases$df)
save_csv(dtdiff_cases, "dtdiff_test_cases",
         "dtdiff (derivative of t-density) for various (x, df)")

# ── 6.  Simulation (with fixed seed) ──────────────────────────────────────
cat("\n[Step 6] Simulation\n")

# Use existing simulation.rds from the R pipeline
# But also generate a NEW stationary simulation with fixed R seed
set.seed(42)
sim_stat <- sim_expt_2d(X, Y, df_true, alpha_true, a=2, b=1, g=0, n_sim=5)
# Save as flat CSV: columns = sim1..sim5, rows = grid points in column-major order
sim_stat_flat <- matrix(NA, nrow = n_grid, ncol = 5)
for (k in 1:5) {
  sim_stat_flat[, k] <- as.vector(sim_stat[,,k])  # column-major
}
colnames(sim_stat_flat) <- paste0("sim", 1:5)
sim_stat_df <- data.frame(flat_index = 1:n_grid, sim_stat_flat)
save_csv(sim_stat_df, "simulation_stationary_seed42",
         "Stationary simulation (a=2,b=1,g=0,df=5,alpha=1, seed=42, 5 sims)")

# Non-stationary simulation with stripes parameters
set.seed(42)
sim_nonstat <- sim_expt_2d_nonstat(X, Y, df_true, alpha_true,
                                    a_matrix, b_matrix, g_matrix, n_sim=10)
sim_nonstat_flat <- matrix(NA, nrow = n_grid, ncol = 10)
for (k in 1:10) {
  sim_nonstat_flat[, k] <- as.vector(sim_nonstat[,,k])
}
colnames(sim_nonstat_flat) <- paste0("sim", 1:10)
sim_nonstat_df <- data.frame(flat_index = 1:n_grid, sim_nonstat_flat)
save_csv(sim_nonstat_df, "simulation_nonstat_seed42",
         "Non-stationary stripes simulation (seed=42, 10 sims)")

# Also save the covariance matrix used for the stationary simulation
cov_matrix_stat <- matrix(sapply(1:n_grid, function(i) {
  sapply(1:n_grid, function(j) {
    cov_fkt_2d(X[i]-X[j], Y[i]-Y[j], alpha_true, 2, 1, 0)
  })
}), nrow = n_grid)
save_csv(as.data.frame(cov_matrix_stat), "covariance_matrix_stationary",
         "100x100 covariance matrix for stationary sim (a=2,b=1,g=0)")

# And the non-stationary covariance matrix
cov_matrix_nonstat <- matrix(sapply(1:n_grid, function(i) {
  sapply(1:n_grid, function(j) {
    cov_fkt_2d_nonstat2(X[i]-X[j], Y[i]-Y[j], alpha_true,
                        a_matrix[i], b_matrix[i], g_matrix[i],
                        a_matrix[j], b_matrix[j], g_matrix[j])
  })
}), nrow = n_grid)
save_csv(as.data.frame(cov_matrix_nonstat), "covariance_matrix_nonstationary",
         "100x100 non-stationary covariance matrix for stripes preset")

# ── 7.  Madogram / Extremal Coefficient Matrix ────────────────────────────
cat("\n[Step 7] Madogram / Extremal Coefficients\n")

# Use the non-stationary simulation for the madogram
a_sim_for_mado <- sim_nonstat  # 10x10x10
# Also store the R n_sim for the functions that need it
n_sim_save <- n_sim  # functions.R uses global n_sim

# Madogram matrix
mado_matrix <- c_extrcoeff_matrix(a_sim_for_mado, madogram = TRUE)
save_csv(as.data.frame(mado_matrix), "madogram_matrix",
         "100x100 madogram (nu_ij) matrix from stripes simulation")

# Extremal coefficient matrix (theta - 1)
ec_matrix <- c_extrcoeff_matrix(a_sim_for_mado, madogram = FALSE)
save_csv(as.data.frame(ec_matrix), "extremal_coefficient_matrix",
         "100x100 extremal coefficient (theta-1) matrix")

# ── 8.  Local Estimation (one grid point) ─────────────────────────────────
cat("\n[Step 8] Local Estimation (single point)\n")

# Run local estimation at a few specific grid points using the non-stat sim
# Use fixed seed for reproducibility of the LHS starting points
a_sim_exp_ns <- sim_nonstat  # global variable expected by functions

locest_results <- matrix(NA, nrow = 5, ncol = 3)
test_points <- c(1, 25, 50, 55, 100)  # R 1-based indices into grid
for (idx in 1:length(test_points)) {
  pt <- test_points[idx]
  x_pt <- X[pt]
  y_pt <- Y[pt]
  cat("  Estimating at grid point", pt, "(x=", x_pt, ", y=", y_pt, ")\n")
  set.seed(42)  # reproducible LHS
  est <- pairwise_density_optim_local(
    a_sim_exp_ns, df_true, alpha_true, x_pt, y_pt,
    abstand = locest_abst, print = FALSE, ensemble = locest_ensemble,
    lower_bounds = c(0.01, 0.01),
    upper_bounds = c(max(a_matrix) + 5, 2 * max(b_matrix))
  )
  locest_results[idx, ] <- est
}
locest_df <- data.frame(
  grid_point = test_points,
  X = X[test_points],
  Y = Y[test_points],
  a_est = locest_results[, 1],
  b_est = locest_results[, 2],
  g_est = locest_results[, 3],
  a_true = a_matrix[test_points],
  b_true = b_matrix[test_points],
  g_true = g_matrix[test_points]
)
save_csv(locest_df, "local_estimates_selected_points",
         "Local estimates at 5 grid points (R 1-based), with true values")

# ── 9.  Full Local Estimation (all grid points) ───────────────────────────
cat("\n[Step 9] Full Local Estimation (all points)\n")

full_locest <- matrix(NA, nrow = n_grid, ncol = 3)
upper_bounds_locest <- c(max(a_matrix) + 5, 2 * max(b_matrix))
for (pt in 1:n_grid) {
  x_pt <- X[pt]
  y_pt <- Y[pt]
  set.seed(42)
  full_locest[pt, ] <- pairwise_density_optim_local(
    a_sim_exp_ns, df_true, alpha_true, x_pt, y_pt,
    abstand = locest_abst, print = FALSE, ensemble = locest_ensemble,
    lower_bounds = c(0.01, 0.01), upper_bounds = upper_bounds_locest
  )
}
full_locest_df <- data.frame(
  flat_index = 1:n_grid,
  a_est = full_locest[, 1],
  b_est = full_locest[, 2],
  g_est = full_locest[, 3]
)
save_csv(full_locest_df, "local_estimates_all",
         "Local estimates for all 100 grid points (raw, before smoothing)")

# ── 10.  Smoothed Local Estimates ─────────────────────────────────────────
cat("\n[Step 10] Smoothing\n")

smoothed <- smooth_local_estimates(full_locest, smoothing_dist)
smoothed_df <- data.frame(
  flat_index = 1:n_grid,
  a_sm = smoothed[, 1],
  b_sm = smoothed[, 2],
  g_sm = smoothed[, 3]
)
save_csv(smoothed_df, "local_estimates_smoothed",
         "Smoothed local estimates (smoothing_dist=1)")

# ── 11.  Ellipse Dissimilarity Matrix (LEC Distance) ─────────────────────
cat("\n[Step 11] Ellipse Dissimilarity\n")

ell_dist <- calc_distance_ellipses(smoothed, res = 21)
save_csv(as.data.frame(ell_dist), "ellipse_dissimilarity_matrix",
         "100x100 LEC ellipse dissimilarity matrix (x100)")

# ── 12.  Hierarchical Clustering (LEC) ───────────────────────────────────
cat("\n[Step 12] LEC Clustering\n")

hc_lec <- clustering(ell_dist)
# Extract merge matrix and heights for comparison
hc_lec_df <- data.frame(
  merge_i = hc_lec$merge[, 1],
  merge_j = hc_lec$merge[, 2],
  height  = hc_lec$height
)
save_csv(hc_lec_df, "hclust_lec_details",
         "LEC hierarchical clustering: merge indices, heights")
save_csv(data.frame(order = hc_lec$order), "hclust_lec_order",
         "LEC hierarchical clustering: leaf order")

# Determine k via threshold method (10th percentile)
q10_lec <- quantile(ell_dist[upper.tri(ell_dist)], 0.10)
k_lec   <- cluster_number_threshold_method(hc_lec, q10_lec)
if (k_lec < 5) k_lec <- 5  # fallback from R code

clusters_lec <- cutree(hc_lec, k_lec)
cluster_lec_df <- data.frame(
  flat_index = 1:n_grid,
  cluster = clusters_lec
)
save_csv(cluster_lec_df, "clusters_lec",
         paste0("LEC cluster assignments (k=", k_lec, ", threshold=q10)"))

save_csv(data.frame(
  k_lec = k_lec,
  q10_threshold = q10_lec
), "lec_clustering_params",
   "LEC clustering parameters: k and q10 threshold")

# ── 13.  Saunders (EDC) Clustering ────────────────────────────────────────
cat("\n[Step 13] EDC Clustering\n")

q30_saunders <- quantile(mado_matrix[upper.tri(mado_matrix)], 0.30)
hc_saunders  <- clustering(mado_matrix)

hc_saunders_df <- data.frame(
  merge_i = hc_saunders$merge[, 1],
  merge_j = hc_saunders$merge[, 2],
  height  = hc_saunders$height
)
save_csv(hc_saunders_df, "hclust_saunders_details",
         "EDC hierarchical clustering: merge indices, heights")
save_csv(data.frame(order = hc_saunders$order), "hclust_saunders_order",
         "EDC hierarchical clustering: leaf order")

k_saunders <- cluster_number_threshold_method(hc_saunders, q30_saunders)
if (k_saunders < 5) k_saunders <- 5

clusters_saunders <- cutree(hc_saunders, k_saunders)
cluster_saunders_df <- data.frame(
  flat_index = 1:n_grid,
  cluster = clusters_saunders
)
save_csv(cluster_saunders_df, "clusters_saunders",
         paste0("Saunders/EDC cluster assignments (k=", k_saunders, ", threshold=q30)"))

save_csv(data.frame(
  k_saunders = k_saunders,
  q30_threshold = q30_saunders
), "saunders_clustering_params",
   "Saunders clustering parameters: k and q30 threshold")

# ── 14.  Log-Likelihood in Clusters (LEC) ────────────────────────────────
cat("\n[Step 14] Log-likelihood in LEC Clusters\n")

llh_lec_results <- data.frame(
  cluster = integer(0),
  a = numeric(0), b = numeric(0), g = numeric(0),
  n_cells = integer(0),
  llh = numeric(0)
)

for (cl in 1:k_lec) {
  which_cl <- which(clusters_lec == cl)
  if (length(which_cl) >= 2) {
    # Use smoothed estimates of the cluster centroid as representative
    # For llh_in_cluster, use the mean of the cluster's smoothed estimates
    mean_est <- colMeans(smoothed[which_cl, , drop = FALSE])
    lh <- llh_in_cluster(
      t(sapply(which_cl, function(j) {a_sim_for_mado[number_grid(j)[1], number_grid(j)[2], ]})),
      df_true, alpha_true,
      X[which_cl], Y[which_cl],
      mean_est,
      max_dist = 4 * ((max(X) - min(X)) / (resolution - 1)),
      average = TRUE
    )
    llh_lec_results <- rbind(llh_lec_results, data.frame(
      cluster = cl, a = mean_est[1], b = mean_est[2], g = mean_est[3],
      n_cells = length(which_cl), llh = lh
    ))
  }
}
save_csv(llh_lec_results, "llh_in_clusters_lec",
         "Log-likelihood per LEC cluster using mean smoothed estimates")

# ── 15.  Summary ──────────────────────────────────────────────────────────
cat("\n[Step 15] Summary\n")

# grid_number / number_grid test cases
coord_tests <- data.frame(
  i = c(1, 5, 10, 3, 7),
  j = c(1, 5, 10, 8, 2)
)
coord_tests$grid_num <- mapply(grid_number, coord_tests$i, coord_tests$j)
save_csv(coord_tests, "grid_number_test_cases",
         "grid_number(i,j) results for verification (R 1-based)")

# number_grid inverse
ncoord_tests <- data.frame(n = c(1, 10, 45, 50, 100))
ng_results <- t(sapply(ncoord_tests$n, number_grid))
ncoord_tests$i <- ng_results[, 1]
ncoord_tests$j <- ng_results[, 2]
save_csv(ncoord_tests, "number_grid_test_cases",
         "number_grid(n) results for verification (R 1-based)")

# koord_num test cases
kn_tests <- data.frame(
  x = c(-5, 0, 5, -2.5, 3.3),
  y = c(5, 0, -5, 2.5, -1.1)
)
kn_tests$grid_num <- mapply(koord_num, kn_tests$x, kn_tests$y)
kn_tests$Xval <- X[kn_tests$grid_num]
kn_tests$Yval <- Y[kn_tests$grid_num]
save_csv(kn_tests, "koord_num_test_cases",
         "koord_num(x,y) nearest-grid-point results (R 1-based index)")

close(manifest)

cat("\n=== Done! Reference data saved to", outdir, "===\n")
cat("Files generated:\n")
system(paste("ls -la", outdir, "| grep .csv"))

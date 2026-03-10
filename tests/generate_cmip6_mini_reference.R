#!/usr/bin/env Rscript
###############################################################################
# Generate reference data for a miniature CMIP6-style pipeline.
#
# Runs a small-scale version of the full Figure 9 workflow in R:
#   monthly precip → STL detrend → annual maxima → GEV → unit Fréchet
#   → EDC madogram → clustering → k
#
# Saves all intermediate results in tests/reference_data/cmip6_mini/
# for verification by Python tests (test_cmip6_rparity.py).
#
# Usage:  cd weatherisk && Rscript tests/generate_cmip6_mini_reference.R
###############################################################################

cat("=== Mini CMIP6 Pipeline Reference Data Generator ===\n")
cat("Working directory:", getwd(), "\n\n")

# ── Load R functions (but NOT parameters — we set our own) ──────────────
source("r_code/functions.R")
library(lhs)

outdir <- "tests/reference_data/cmip6_mini"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save_csv <- function(obj, name, description) {
  path <- file.path(outdir, paste0(name, ".csv"))
  write.csv(obj, path, row.names = FALSE)
  cat("  Saved:", name, ".csv  —  ", description, "\n")
}

# ── Override globals needed by functions.R ──────────────────────────────
resolution <- 5
n_years <- 20

x_ax  <- seq(-2, 2, length.out = resolution)
y_ax  <- -seq(-2, 2, length.out = resolution)
X     <- rep(1, length(y_ax)) %*% t(x_ax)
Y     <- y_ax %*% t(rep(1, length(x_ax)))

nrow_grid <- length(y_ax)
ncol_grid <- length(x_ax)
n_grid    <- nrow_grid * ncol_grid
exp.grid  <- expand.grid(y = y_ax, x = x_ax)

df_true      <- 5
alpha_true   <- 1
n_sim        <- n_years   # c_extrcoeff_matrix uses global n_sim

cat("Grid:", nrow_grid, "x", ncol_grid, "=", n_grid, "cells\n")
cat("Years:", n_years, "(", n_years * 12, "months)\n\n")

# ── 1. Generate synthetic monthly precipitation ───────────────────────
cat("[Step 1] Generating synthetic monthly precipitation\n")

set.seed(42)
n_months <- n_years * 12
pr <- array(0, dim = c(n_months, nrow_grid, ncol_grid))

for (i in 1:nrow_grid) {
  for (j in 1:ncol_grid) {
    t_m <- 1:n_months
    amp   <- 3 + 2 * (j / ncol_grid)
    trend <- 0.001 * (1 + i / nrow_grid) * t_m
    seasonal <- amp * sin(2 * pi * t_m / 12)
    noise <- rnorm(n_months, 0, 1.5)
    pr[, i, j] <- pmax(15 + seasonal + trend + noise, 0.1)
  }
}

# Save as (n_months, n_cells) with cells in column-major R order
pr_flat <- matrix(pr, nrow = n_months, ncol = n_grid)
save_csv(as.data.frame(pr_flat), "monthly_pr",
         paste0("Raw monthly precip (", n_months, " x ", n_grid, ") col-major"))

# Grid coordinates (column-major R order)
grid_df <- data.frame(
  flat_index = 1:n_grid,
  X = as.vector(X),
  Y = as.vector(Y),
  row_i = rep(1:nrow_grid, ncol_grid),   # R 1-based
  col_j = rep(1:ncol_grid, each = nrow_grid)
)
save_csv(grid_df, "grid_coordinates",
         "Grid coordinates in column-major R order")

# Save scalar parameters
save_csv(data.frame(
  parameter = c("resolution", "n_years", "n_months", "nrow_grid", "ncol_grid",
                 "n_grid", "df_true", "alpha_true"),
  value = c(resolution, n_years, n_months, nrow_grid, ncol_grid,
            n_grid, df_true, alpha_true)
), "scalar_params", "Scalar parameters")

# ── 2. STL detrending ────────────────────────────────────────────────
cat("\n[Step 2] STL detrending\n")

detrended  <- pr
stl_trends <- array(0, dim = c(n_months, nrow_grid, ncol_grid))

for (i in 1:nrow_grid) {
  for (j in 1:ncol_grid) {
    ts_data <- ts(pr[, i, j], frequency = 12)
    stl_fit <- stl(ts_data, s.window = "periodic", robust = TRUE)
    stl_trends[, i, j] <- stl_fit$time.series[, "trend"]
    detrended[, i, j]  <- pr[, i, j] - stl_trends[, i, j]
  }
}

save_csv(as.data.frame(matrix(detrended,  nrow = n_months, ncol = n_grid)),
         "detrended", "STL detrended monthly data")
save_csv(as.data.frame(matrix(stl_trends, nrow = n_months, ncol = n_grid)),
         "stl_trends", "STL trend component")

# ── 3. Annual maxima ─────────────────────────────────────────────────
cat("\n[Step 3] Annual maxima of detrended data\n")

am <- array(0, dim = c(n_years, nrow_grid, ncol_grid))
for (yr in 1:n_years) {
  idx <- ((yr - 1) * 12 + 1):(yr * 12)
  am[yr, , ] <- apply(detrended[idx, , , drop = FALSE], c(2, 3), max)
}

save_csv(as.data.frame(matrix(am, nrow = n_years, ncol = n_grid)),
         "annual_maxima", "Annual maxima of detrended data")

# ── 4. GEV fit + Fréchet transform ──────────────────────────────────
cat("\n[Step 4] GEV fitting and Fréchet transform\n")

frechet    <- am
gev_params <- data.frame(cell = 1:n_grid, loc = NA, scale = NA, shape = NA)

for (s in 1:n_grid) {
  i <- ((s - 1) %%  nrow_grid) + 1
  j <- ((s - 1) %/% nrow_grid) + 1
  y_data <- am[, i, j]

  # Simple MLE GEV using optim (avoid external package dependency)
  gev_nll <- function(par) {
    mu <- par[1]; sigma <- par[2]; xi <- par[3]
    if (sigma <= 0) return(1e10)
    z <- (y_data - mu) / sigma
    if (abs(xi) < 1e-8) {
      return(sum(log(sigma) + z + exp(-z)))
    }
    w <- 1 + xi * z
    if (any(w <= 0)) return(1e10)
    return(sum(log(sigma) + (1 + 1/xi) * log(w) + w^(-1/xi)))
  }

  # Starting values from method of moments
  mu0    <- mean(y_data)
  sig0   <- sd(y_data) * sqrt(6) / pi
  init   <- c(mu0 - 0.5772 * sig0, sig0, 0.1)
  result <- optim(init, gev_nll, method = "Nelder-Mead",
                  control = list(maxit = 5000))

  mu  <- result$par[1]
  sig <- result$par[2]
  xi  <- result$par[3]

  gev_params$loc[s]   <- mu
  gev_params$scale[s] <- sig
  gev_params$shape[s] <- xi

  # Transform to unit Fréchet
  if (abs(xi) < 1e-8) {
    frechet[, i, j] <- exp((y_data - mu) / sig)
  } else {
    w <- 1 + xi * (y_data - mu) / sig
    w <- pmax(w, 1e-10)
    frechet[, i, j] <- w^(1 / xi)
  }
}

save_csv(gev_params, "gev_params", "GEV fit parameters per cell")
save_csv(as.data.frame(matrix(frechet, nrow = n_years, ncol = n_grid)),
         "frechet", "Fréchet-transformed annual maxima")

# ── 5. EDC madogram + clustering ─────────────────────────────────────
cat("\n[Step 5] EDC madogram and clustering\n")

# c_extrcoeff_matrix expects (nrow, ncol, n_sim) — transpose from (n_years, nrow, ncol)
a_sim_local <- aperm(frechet, c(2, 3, 1))
v_matrix <- c_extrcoeff_matrix(a_sim_local, madogram = TRUE)

save_csv(as.data.frame(v_matrix), "madogram_matrix",
         paste0(n_grid, "x", n_grid, " madogram matrix"))

# Cluster
v_tri      <- v_matrix[upper.tri(v_matrix)]
q30        <- quantile(v_tri, 0.3)
hc_edc     <- hclust(as.dist(v_matrix), method = "average")
k_edc      <- cluster_number_threshold_method(hc_edc, q30)
if (k_edc < 2) k_edc <- 2
clust_edc  <- cutree(hc_edc, k_edc)

save_csv(data.frame(k_edc = k_edc, q30 = q30),
         "edc_clustering_params", "EDC k and 30th-percentile threshold")
save_csv(data.frame(flat_index = 1:n_grid, cluster = clust_edc),
         "clusters_edc", paste0("EDC cluster assignments (k=", k_edc, ")"))

cat("  k_EDC =", k_edc, "  q30 =", round(q30, 5), "\n")

# ── 6. Local MLE estimation (all grid points) ───────────────────────
cat("\n[Step 6] Local MLE estimation\n")

locest_ensemble <- 5
locest_abst     <- 2   # neighbor radius in grid-point units
smoothing_dist  <- 2

full_locest <- matrix(NA, nrow = n_grid, ncol = 3)
for (pt in 1:n_grid) {
  x_pt <- X[pt]
  y_pt <- Y[pt]
  set.seed(42)
  full_locest[pt, ] <- pairwise_density_optim_local(
    a_sim_local, df_true, alpha_true, x_pt, y_pt,
    abstand = locest_abst, print = FALSE, ensemble = locest_ensemble,
    lower_bounds = c(0.01, 0.01), upper_bounds = c(15, 15)
  )
  cat("  Cell", pt, "/", n_grid,
      ": a=", round(full_locest[pt, 1], 3),
      " b=", round(full_locest[pt, 2], 3),
      " g=", round(full_locest[pt, 3], 3), "\n")
}

save_csv(data.frame(
  flat_index = 1:n_grid,
  a_est = full_locest[, 1],
  b_est = full_locest[, 2],
  g_est = full_locest[, 3]
), "local_estimates_all", "Local MLE estimates for all grid points (raw)")

# ── 7. Spatial smoothing ────────────────────────────────────────────
cat("\n[Step 7] Spatial smoothing (smoothing_dist =", smoothing_dist, ")\n")

smoothed <- smooth_local_estimates(full_locest, smoothing_dist)

save_csv(data.frame(
  flat_index = 1:n_grid,
  a_sm = smoothed[, 1],
  b_sm = smoothed[, 2],
  g_sm = smoothed[, 3]
), "local_estimates_smoothed", "Smoothed local estimates")

# ── 8. Ellipse dissimilarity (LEC distance) ────────────────────────
cat("\n[Step 8] Ellipse dissimilarity\n")

ell_dist <- calc_distance_ellipses(smoothed, res = 21)

save_csv(as.data.frame(ell_dist), "ellipse_dissimilarity_matrix",
         paste0(n_grid, "x", n_grid, " LEC ellipse dissimilarity (x100)"))

# ── 9. LEC clustering ──────────────────────────────────────────────
cat("\n[Step 9] LEC clustering\n")

hc_lec   <- clustering(ell_dist)
q30_lec  <- quantile(ell_dist[upper.tri(ell_dist)], 0.30)
k_lec    <- cluster_number_threshold_method(hc_lec, q30_lec)
if (k_lec < 2) k_lec <- 2
clust_lec <- cutree(hc_lec, k_lec)

save_csv(data.frame(k_lec = k_lec, q30_lec = q30_lec),
         "lec_clustering_params", "LEC k and 30th-percentile threshold")
save_csv(data.frame(flat_index = 1:n_grid, cluster = clust_lec),
         "clusters_lec", paste0("LEC cluster assignments (k=", k_lec, ")"))

cat("  k_LEC =", k_lec, "  q30 =", round(q30_lec, 5), "\n")

# Save LEC-related scalar params
save_csv(data.frame(
  parameter = c("locest_ensemble", "locest_abst", "smoothing_dist"),
  value = c(locest_ensemble, locest_abst, smoothing_dist)
), "lec_scalar_params", "LEC pipeline parameters")

# ── Summary ──────────────────────────────────────────────────────────
cat("\n=== Mini CMIP6 Reference Data Complete ===\n")
cat("Output directory:", outdir, "\n")
cat("Files: monthly_pr, grid_coordinates, scalar_params,\n")
cat("       detrended, stl_trends, annual_maxima,\n")
cat("       gev_params, frechet, madogram_matrix,\n")
cat("       edc_clustering_params, clusters_edc,\n")
cat("       local_estimates_all, local_estimates_smoothed,\n")
cat("       ellipse_dissimilarity_matrix,\n")
cat("       lec_clustering_params, clusters_lec, lec_scalar_params\n")

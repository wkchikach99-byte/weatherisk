#!/usr/bin/env Rscript

suppressMessages({
  library(lhs)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/replay_fig9_benchmark_in_r.R <input_dir> <output_json>")
}

input_dir <- args[1]
output_json <- args[2]

source("r_code/functions.R")

calc_distance_ellipses_corrected <- function(local_estimates_matrix, res = 21) {
  xs <- rep(seq(-1, 1, length.out = res), res)
  ys <- rep(seq(-1, 1, length.out = res), each = res)
  wh <- which(((xs^2 + ys^2) <= res^2) & (ys > 0) | (xs > 0))
  xs <- xs[wh]
  ys <- ys[wh]

  n <- nrow(local_estimates_matrix)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      mx <- max(
        local_estimates_matrix[i, 1] + local_estimates_matrix[i, 2],
        local_estimates_matrix[j, 1] + local_estimates_matrix[j, 2]
      )
      ell1 <- 1 * (cov_fkt_2d(
        xs,
        ys,
        alpha = 1,
        a = local_estimates_matrix[i, 1] / mx,
        b = local_estimates_matrix[i, 2] / mx,
        g = local_estimates_matrix[i, 3]
      ) > exp(-1))
      ell2 <- 1 * (cov_fkt_2d(
        xs,
        ys,
        alpha = 1,
        a = local_estimates_matrix[j, 1] / mx,
        b = local_estimates_matrix[j, 2] / mx,
        g = local_estimates_matrix[j, 3]
      ) > exp(-1))
      if (max(sum(ell1), sum(ell2)) == 0) {
        dist_matrix[i, j] <- 1
      } else {
        dist_matrix[i, j] <- 1 - (sum(ell1 * ell2) + 1 / 2) / (sum(ell1 + ell2 - ell1 * ell2) + 1 / 2)
      }
    }
  }
  100 * (dist_matrix + t(dist_matrix))
}

read_scalar <- function(df, name) {
  as.numeric(df$value[df$parameter == name][1])
}

scalar_params <- read.csv(file.path(input_dir, "scalar_params.csv"), stringsAsFactors = FALSE)
monthly_pr <- as.matrix(read.csv(file.path(input_dir, "monthly_pr.csv"), check.names = FALSE))

resolution <- as.integer(read_scalar(scalar_params, "resolution"))
n_years <- as.integer(read_scalar(scalar_params, "n_years"))
n_months <- as.integer(read_scalar(scalar_params, "n_months"))
nrow_grid <- resolution
ncol_grid <- resolution
n_grid <- nrow_grid * ncol_grid

df_true <- read_scalar(scalar_params, "df")
alpha_true <- read_scalar(scalar_params, "alpha")
locest_abst <- read_scalar(scalar_params, "neighbor_radius")
smoothing_dist <- read_scalar(scalar_params, "smoothing_radius")
locest_ensemble <- as.integer(read_scalar(scalar_params, "mle_ensemble"))

n_sim <- n_years
x_ax <- seq(0, ncol_grid - 1, by = 1)
y_ax <- seq(0, nrow_grid - 1, by = 1)
X <- rep(1, length(y_ax)) %*% t(x_ax)
Y <- y_ax %*% t(rep(1, length(x_ax)))
exp.grid <- expand.grid(y = y_ax, x = x_ax)

pr <- array(monthly_pr, dim = c(n_months, nrow_grid, ncol_grid))

detrended <- pr
stl_trends <- array(0, dim = c(n_months, nrow_grid, ncol_grid))
for (i in 1:nrow_grid) {
  for (j in 1:ncol_grid) {
    ts_data <- ts(pr[, i, j], frequency = 12)
    stl_fit <- stl(ts_data, s.window = "periodic", robust = TRUE)
    stl_trends[, i, j] <- stl_fit$time.series[, "trend"]
    detrended[, i, j] <- pr[, i, j] - stl_trends[, i, j]
  }
}

am <- array(0, dim = c(n_years, nrow_grid, ncol_grid))
for (yr in 1:n_years) {
  idx <- ((yr - 1) * 12 + 1):(yr * 12)
  am[yr, , ] <- apply(detrended[idx, , , drop = FALSE], c(2, 3), max)
}

frechet <- am
gev_params <- matrix(NA_real_, nrow = n_grid, ncol = 3)
colnames(gev_params) <- c("loc", "scale", "shape")

for (s in 1:n_grid) {
  i <- ((s - 1) %% nrow_grid) + 1
  j <- ((s - 1) %/% nrow_grid) + 1
  y_data <- am[, i, j]

  gev_nll <- function(par) {
    mu <- par[1]
    sigma <- par[2]
    xi <- par[3]
    if (sigma <= 0) return(1e10)
    z <- (y_data - mu) / sigma
    if (abs(xi) < 1e-8) {
      return(sum(log(sigma) + z + exp(-z)))
    }
    w <- 1 + xi * z
    if (any(w <= 0)) return(1e10)
    sum(log(sigma) + (1 + 1 / xi) * log(w) + w^(-1 / xi))
  }

  mu0 <- mean(y_data)
  sig0 <- sd(y_data) * sqrt(6) / pi
  init <- c(mu0 - 0.5772 * sig0, sig0, 0.1)
  fit <- optim(init, gev_nll, method = "Nelder-Mead", control = list(maxit = 5000))

  mu <- fit$par[1]
  sig <- fit$par[2]
  xi <- fit$par[3]
  gev_params[s, ] <- c(mu, sig, xi)

  if (abs(xi) < 1e-8) {
    frechet[, i, j] <- exp((y_data - mu) / sig)
  } else {
    w <- 1 + xi * (y_data - mu) / sig
    w <- pmax(w, 1e-10)
    frechet[, i, j] <- w^(1 / xi)
  }
}

a_sim_local <- aperm(frechet, c(2, 3, 1))
v_matrix <- c_extrcoeff_matrix(a_sim_local, madogram = TRUE)
v_tri <- v_matrix[upper.tri(v_matrix)]
q30_edc <- quantile(v_tri, 0.30)
hc_edc <- hclust(as.dist(v_matrix), method = "average")
k_edc <- cluster_number_threshold_method(hc_edc, q30_edc)
if (k_edc < 2) k_edc <- 2

full_locest <- matrix(NA_real_, nrow = n_grid, ncol = 3)
for (pt in 1:n_grid) {
  x_pt <- X[pt]
  y_pt <- Y[pt]
  set.seed(42)
  full_locest[pt, ] <- pairwise_density_optim_local(
    a_sim_local,
    df_true,
    alpha_true,
    x_pt,
    y_pt,
    abstand = locest_abst,
    print = FALSE,
    ensemble = locest_ensemble,
    lower_bounds = c(0.01, 0.01),
    upper_bounds = c(15, 15)
  )
}

smoothed <- smooth_local_estimates(full_locest, smoothing_dist)
ell_dist_raw <- calc_distance_ellipses(smoothed, res = 21)
ell_dist_corrected <- calc_distance_ellipses_corrected(smoothed, res = 21)

hc_lec_raw <- clustering(ell_dist_raw)
q30_lec_raw <- quantile(ell_dist_raw[upper.tri(ell_dist_raw)], 0.30)
k_lec_raw <- cluster_number_threshold_method(hc_lec_raw, q30_lec_raw)
if (k_lec_raw < 2) k_lec_raw <- 2

hc_lec_corrected <- clustering(ell_dist_corrected)
q30_lec_corrected <- quantile(ell_dist_corrected[upper.tri(ell_dist_corrected)], 0.30)
k_lec_corrected <- cluster_number_threshold_method(hc_lec_corrected, q30_lec_corrected)
if (k_lec_corrected < 2) k_lec_corrected <- 2

result_lines <- c(
  "{",
  sprintf('  "resolution": %d,', resolution),
  sprintf('  "n_years": %d,', n_years),
  sprintf('  "k_edc": %d,', k_edc),
  sprintf('  "k_lec_raw": %d,', k_lec_raw),
  sprintf('  "k_lec_corrected": %d,', k_lec_corrected),
  sprintf('  "q30_edc": %.17g,', q30_edc),
  sprintf('  "q30_lec_raw": %.17g,', q30_lec_raw),
  sprintf('  "q30_lec_corrected": %.17g', q30_lec_corrected),
  "}"
)
writeLines(result_lines, output_json)
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

args <- commandArgs(trailingOnly = TRUE)
optimizer_only <- (
  any(args %in% c("optimizer", "local-optimizer", "--only=optimizer")) ||
  any(head(args, -1) == "--only" & tail(args, -1) == "optimizer")
)
optimizer_min <- any(args %in% c("optimizer-min", "local-optimizer-min", "--only=optimizer-min"))
optimizer_micro <- any(args %in% c("optimizer-micro", "local-optimizer-micro", "--only=optimizer-micro"))

# ── 0.  Source R code and parameters ────────────────────────────────────────
library(lhs)                      # needed by pairwise_density_optim_local
source("r_code/functions.R")      # all math functions
source("r_code/parameters.R")     # stripes preset (resolution=10, etc.)

configure_small_reference_preset <- function() {
  resolution <<- 10
  df_true <<- 5
  alpha_true <<- 1
  n_sim <<- 10

  x_ax <<- seq(-5, 5, length.out = resolution)
  y_ax <<- -seq(-5, 5, length.out = resolution)

  X <<- rep(1, length(y_ax)) %*% t(x_ax)
  Y <<- y_ax %*% t(rep(1, length(x_ax)))

  nrow_grid <<- length(y_ax)
  ncol_grid <<- length(x_ax)
  n_grid <<- nrow_grid * ncol_grid
  exp.grid <<- expand.grid(y = y_ax, x = x_ax)

  b_end <<- 5
  a_matrix <<- 0 * X + 2
  b_matrix <<- (X + 5) / 10 * b_end
  g_matrix <<- 0 * X

  locest_ensemble <<- 5
  locest_abst <<- 4
  smoothing_dist <<- 2
}

if (optimizer_only || optimizer_min) {
  configure_small_reference_preset()
}

# Always configure the small preset for the full reference run too
configure_small_reference_preset()

outdir <- "tests/reference_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(outdir, "MANIFEST.txt")

# Write manifest header
manifest_mode <- if ((optimizer_only || optimizer_min || optimizer_micro) && file.exists(manifest_path)) "a" else "w"
manifest <- file(manifest_path, open = manifest_mode)
if (manifest_mode == "w") {
  writeLines("# R Reference Data Manifest", manifest)
  writeLines(paste0("# Generated: ", Sys.time()), manifest)
  writeLines(paste0("# R version: ", R.version.string), manifest)
  writeLines("", manifest)
} else {
  writeLines("", manifest)
  writeLines(paste0("# Optimizer-only refresh: ", Sys.time()), manifest)
}

save_csv <- function(obj, name, description) {
  path <- file.path(outdir, paste0(name, ".csv"))
  # write.csv + options(digits=22) still caps at ~15 significant digits.
  # Use a custom writer with sprintf("%.17e") for full double precision.
  df <- as.data.frame(obj)
  header <- paste(paste0('"', names(df), '"'), collapse = ",")
  lines <- apply(df, 1, function(row) {
    vals <- sapply(seq_along(row), function(i) {
      v <- row[[i]]
      if (is.na(v)) {
        "NA"
      } else if (is.integer(v) || (is.numeric(v) && v == round(v) && abs(v) < 1e15)) {
        format(v, scientific = FALSE)
      } else if (is.numeric(v)) {
        sprintf("%.17e", as.numeric(v))
      } else {
        as.character(v)
      }
    })
    paste(vals, collapse = ",")
  })
  writeLines(c(header, lines), path)
  writeLines(paste0(name, ".csv  —  ", description), manifest)
  cat("  Saved:", name, "\n")
}

build_local_optimizer_fixture <- function(a_sim, x_pt, y_pt, abstand) {
  xy_pos <- which.min(dist_x(x_pt, X)^2 + dist_y(y_pt, Y)^2)[1]
  xs <- rep(-abstand:abstand, (2 * abstand + 1))
  ys <- rep(-abstand:abstand, each = (2 * abstand + 1))
  which.dist <- which(sqrt(xs * xs + ys * ys) <= abstand & (xs * xs + ys * ys) > 0)
  xs <- xs[which.dist]
  ys <- ys[which.dist]

  sel.x <- number_grid(xy_pos)[2] + xs
  sel.y <- number_grid(xy_pos)[1] + ys
  in.bounds <- which(sel.y > 0 & sel.y <= nrow_grid & sel.x > 0 & sel.x <= ncol_grid)
  sel.x <- sel.x[in.bounds]
  sel.y <- sel.y[in.bounds]
  sel.grid <- sapply(seq_along(in.bounds), function(i) { grid_number(sel.y[i], sel.x[i]) })

  n_sel <- length(sel.x)
  zlist1 <- as.vector(sapply(seq_len(n_sel), function(i) {
    a_sim[number_grid(sel.grid[i])[1], number_grid(sel.grid[i])[2], ]
  }))
  zlist2 <- rep(a_sim[number_grid(xy_pos)[1], number_grid(xy_pos)[2], ], n_sel)
  Xlist <- rep(X[sel.grid] - X[xy_pos], each = n_sim)
  Ylist <- rep(Y[sel.grid] - Y[xy_pos], each = n_sim)

  data.frame(
    grid_point = xy_pos,
    X = x_pt,
    Y = y_pt,
    obs_index = seq_along(zlist1),
    zi = zlist1,
    zj = zlist2,
    xl = Xlist,
    yl = Ylist
  )
}

# Run optimizer on pre-built arrays, matching the exact R code used by
# Python's _optimize_local_mle_via_r subprocess bridge.
run_optimizer_on_arrays <- function(zi, zj, xl, yl, df, alpha,
                                   ensemble, seed, max_boundary_retries,
                                   lower_bounds, upper_bounds) {
  llh <- function(par) {
    sum(pairwise_density_summand(zi, zj, xl, yl, df, alpha, par[1], par[2], par[3]))
  }

  parameters_lower_bound <- c(lower_bounds[1], lower_bounds[2], -pi / 2)
  parameters_upper_bound <- c(upper_bounds[1], upper_bounds[2], pi / 2)
  parscale <- (parameters_upper_bound - parameters_lower_bound) / 100

  set.seed(seed)
  starts <- matrix(
    parameters_lower_bound,
    ensemble + max_boundary_retries,
    length(parameters_lower_bound),
    byrow = TRUE
  ) + maximinLHS(
    ensemble + max_boundary_retries,
    length(parameters_lower_bound)
  ) %*% diag(parameters_upper_bound - parameters_lower_bound)

  num_calc <- 1L
  num_more_calc_done <- 0L
  fit <- NULL
  while (num_calc <= ensemble) {
    fit_ <- optim(
      starts[num_calc + num_more_calc_done, ],
      fn = llh, method = "L-BFGS-B",
      lower = parameters_lower_bound, upper = parameters_upper_bound,
      control = list(fnscale = -1, parscale = parscale, maxit = 10000)
    )
    if (abs(fit_$par[3]) == pi / 2) {
      fit_ <- optim(
        c(fit_$par[1], fit_$par[2], -fit_$par[3]),
        fn = llh, method = "L-BFGS-B",
        lower = parameters_lower_bound, upper = parameters_upper_bound,
        control = list(fnscale = -1, parscale = parscale, maxit = 10000)
      )
    }
    if (is.null(fit) || (-fit_$value < -fit$value)) {
      fit <- fit_
    }
    if (min(abs(fit_$par - parameters_lower_bound)) < 0.01 ||
        min(abs(fit_$par - parameters_upper_bound)) < 0.01) {
      if (num_more_calc_done < max_boundary_retries) {
        num_more_calc_done <- num_more_calc_done + 1L
        next
      }
    }
    num_calc <- num_calc + 1L
  }

  return(as.vector(fit$par))
}

write_optimizer_fixtures <- function(test_points, inputs_name, outputs_name, point_label) {
  cat("\n[Fast Path] ", point_label, "\n", sep = "")

  set.seed(42)
  locest_lhs <- maximinLHS(locest_ensemble + 5, 3)
  colnames(locest_lhs) <- c("u1", "u2", "u3")
  save_csv(as.data.frame(locest_lhs), "maximin_lhs_10x3_seed42",
           "R maximinLHS(10,3) start matrix used by local estimation with set.seed(42)")

  set.seed(42)
  a_sim_exp_ns <- sim_expt_2d_nonstat(X, Y, df_true, alpha_true,
                                      a_matrix, b_matrix, g_matrix, n_sim = n_sim)

  optimizer_input_rows <- vector("list", length(test_points))
  upper_a <- max(a_matrix) + 5
  upper_b <- 2 * max(b_matrix)

  # Build inputs and save as both CSV and binary
  for (idx in seq_along(test_points)) {
    pt <- test_points[idx]
    x_pt <- X[pt]
    y_pt <- Y[pt]
    optimizer_input_rows[[idx]] <- build_local_optimizer_fixture(a_sim_exp_ns, x_pt, y_pt, locest_abst)
  }
  all_inputs <- do.call(rbind, optimizer_input_rows)
  save_csv(all_inputs, inputs_name,
           paste0("Lossless local-optimizer pair arrays (zi,zj,xl,yl) for ", point_label))

  # Save as binary for parse-exact transfer
  bin_path <- file.path(outdir, paste0(inputs_name, ".bin"))
  con <- file(bin_path, "wb")
  writeBin(as.integer(nrow(all_inputs)), con, size = 4, endian = "little")
  writeBin(as.double(all_inputs$grid_point), con, size = 8, endian = "little")
  writeBin(as.double(all_inputs$zi), con, size = 8, endian = "little")
  writeBin(as.double(all_inputs$zj), con, size = 8, endian = "little")
  writeBin(as.double(all_inputs$xl), con, size = 8, endian = "little")
  writeBin(as.double(all_inputs$yl), con, size = 8, endian = "little")
  close(con)
  cat("  Saved:", inputs_name, ".bin (binary)\n")

  # Run optimizer on the in-memory arrays (same binary-exact values)
  locest_results <- matrix(NA_real_, nrow = length(test_points), ncol = 3)
  for (idx in seq_along(test_points)) {
    pt <- test_points[idx]
    inp <- all_inputs[all_inputs$grid_point == pt, ]
    locest_results[idx, ] <- run_optimizer_on_arrays(
      inp$zi, inp$zj, inp$xl, inp$yl,
      df_true, alpha_true,
      ensemble = locest_ensemble, seed = 42, max_boundary_retries = 5,
      lower_bounds = c(0.01, 0.01), upper_bounds = c(upper_a, upper_b)
    )
  }

  save_csv(data.frame(
    grid_point = test_points,
    X = X[test_points],
    Y = Y[test_points],
    lower_a = rep(0.01, length(test_points)),
    lower_b = rep(0.01, length(test_points)),
    lower_g = rep(-pi / 2, length(test_points)),
    upper_a = rep(upper_a, length(test_points)),
    upper_b = rep(upper_b, length(test_points)),
    upper_g = rep(pi / 2, length(test_points)),
    ensemble = rep(locest_ensemble, length(test_points)),
    abstand = rep(locest_abst, length(test_points)),
    a_est = locest_results[, 1],
    b_est = locest_results[, 2],
    g_est = locest_results[, 3]
  ), outputs_name,
  paste0("Exact R-selected optimizer outputs and bounds for ", point_label))
}

write_optimizer_micro_fixtures <- function() {
  cat("\n[Fast Path] Optimizer micro fixtures\n")

  micro_cases <- list(
    list(
      case_id = "interior_scaled",
      start = c(4.8, 8.7, 1.1),
      lower = c(0.01, 0.01, -pi / 2),
      upper = c(7.0, 10.0, pi / 2),
      fn = function(par) {
        ((par[1] - 1.25) / 0.35)^2 + ((par[2] - 3.4) / 1.8)^2 + ((par[3] + 0.42) / 0.12)^2
      }
    ),
    list(
      case_id = "upper_boundary_b",
      start = c(6.4, 0.4, -1.1),
      lower = c(0.01, 0.01, -pi / 2),
      upper = c(7.0, 10.0, pi / 2),
      fn = function(par) {
        (par[1] - 0.45)^2 + (par[2] - 10.8)^2 + 0.2 * (par[3] + 0.3)^2
      }
    ),
    list(
      case_id = "tilted_valley",
      start = c(5.9, 9.1, 0.9),
      lower = c(0.01, 0.01, -pi / 2),
      upper = c(7.0, 10.0, pi / 2),
      fn = function(par) {
        dx <- par[1] - 0.55
        dy <- par[2] - 6.0
        dz <- par[3] + 0.65
        dx^2 + 2.0 * dy^2 + 0.5 * dz^2 + 0.6 * dx * dy
      }
    )
  )

  rows <- lapply(micro_cases, function(case) {
    fit <- optim(
      case$start,
      fn = case$fn,
      method = c("L-BFGS-B"),
      lower = case$lower,
      upper = case$upper,
      control = list(
        parscale = (case$upper - case$lower) / 100,
        maxit = 10000
      )
    )

    data.frame(
      case_id = case$case_id,
      start_a = case$start[1],
      start_b = case$start[2],
      start_g = case$start[3],
      lower_a = case$lower[1],
      lower_b = case$lower[2],
      lower_g = case$lower[3],
      upper_a = case$upper[1],
      upper_b = case$upper[2],
      upper_g = case$upper[3],
      opt_a = fit$par[1],
      opt_b = fit$par[2],
      opt_g = fit$par[3],
      opt_value = fit$value
    )
  })

  save_csv(do.call(rbind, rows), "optim_lbfgsb_micro_cases",
           "R optim(L-BFGS-B) results on toy bounded objectives with parscale=(upper-lower)/100")
}

if (optimizer_min) {
  write_optimizer_fixtures(c(1), "local_optimizer_min_inputs", "local_optimizer_min_outputs", "the single selected local-estimation cell (grid point 1)")
  close(manifest)
  quit(save = "no", status = 0)
}

if (optimizer_micro) {
  write_optimizer_micro_fixtures()
  close(manifest)
  quit(save = "no", status = 0)
}

if (optimizer_only) {
  write_optimizer_fixtures(c(1, 25, 50, 55, 100), "local_optimizer_selected_inputs", "local_optimizer_selected_outputs", "the 5 selected local-estimation cells")
  close(manifest)
  quit(save = "no", status = 0)
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

# Also save as binary (raw double bytes) to avoid parser differences between R and Python
sim_nonstat_bin_path <- file.path(outdir, "simulation_nonstat_seed42.bin")
con <- file(sim_nonstat_bin_path, "wb")
# Layout: n_grid doubles for sim1, then n_grid for sim2, ..., n_grid for sim10
for (k in 1:10) {
  writeBin(sim_nonstat_flat[, k], con, size = 8, endian = "little")
}
close(con)
cat("  Saved: simulation_nonstat_seed42.bin (binary)\n")

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
# Override global n_sim (250 from parameters.R) to match our 10-sim fixture
n_sim_save <- n_sim
n_sim <- 10

# Madogram matrix
mado_matrix <- c_extrcoeff_matrix(a_sim_for_mado, madogram = TRUE)
save_csv(as.data.frame(mado_matrix), "madogram_matrix",
         "100x100 madogram (nu_ij) matrix from stripes simulation")
# Read back from CSV so downstream hclust uses exactly the serialized values
mado_matrix <- as.matrix(read.csv(file.path(outdir, "madogram_matrix.csv")))

# Extremal coefficient matrix (theta - 1)
ec_matrix <- c_extrcoeff_matrix(a_sim_for_mado, madogram = FALSE)
save_csv(as.data.frame(ec_matrix), "extremal_coefficient_matrix",
         "100x100 extremal coefficient (theta-1) matrix")

# ── 8.  Local Estimation (one grid point) ─────────────────────────────────
cat("\n[Step 8] Local Estimation (single point)\n")

set.seed(42)
locest_lhs <- maximinLHS(locest_ensemble + 5, 3)
colnames(locest_lhs) <- c("u1", "u2", "u3")
save_csv(as.data.frame(locest_lhs), "maximin_lhs_10x3_seed42",
         "R maximinLHS(10,3) start matrix used by local estimation with set.seed(42)")

# Run local estimation at a few specific grid points using the non-stat sim
a_sim_exp_ns <- sim_nonstat  # global variable expected by functions

test_points <- c(1, 25, 50, 55, 100)  # R 1-based indices into grid
upper_a <- max(a_matrix) + 5
upper_b <- 2 * max(b_matrix)

# Build inputs and save as both CSV and binary
optimizer_input_rows <- list()
for (idx in seq_along(test_points)) {
  pt <- test_points[idx]
  x_pt <- X[pt]
  y_pt <- Y[pt]
  cat("  Building inputs for grid point", pt, "(x=", x_pt, ", y=", y_pt, ")\n")
  optimizer_input_rows[[idx]] <- build_local_optimizer_fixture(a_sim_exp_ns, x_pt, y_pt, locest_abst)
}
all_inputs <- do.call(rbind, optimizer_input_rows)
save_csv(all_inputs, "local_optimizer_selected_inputs",
         "Lossless local-optimizer pair arrays (zi,zj,xl,yl) for the 5 selected grid points")

# Save as binary for parse-exact transfer
bin_path <- file.path(outdir, "local_optimizer_selected_inputs.bin")
con <- file(bin_path, "wb")
writeBin(as.integer(nrow(all_inputs)), con, size = 4, endian = "little")
writeBin(as.double(all_inputs$grid_point), con, size = 8, endian = "little")
writeBin(as.double(all_inputs$zi), con, size = 8, endian = "little")
writeBin(as.double(all_inputs$zj), con, size = 8, endian = "little")
writeBin(as.double(all_inputs$xl), con, size = 8, endian = "little")
writeBin(as.double(all_inputs$yl), con, size = 8, endian = "little")
close(con)
cat("  Saved: local_optimizer_selected_inputs.bin (binary)\n")

# Run optimizer on the in-memory arrays (binary-exact values)
locest_results <- matrix(NA, nrow = 5, ncol = 3)
for (idx in seq_along(test_points)) {
  pt <- test_points[idx]
  cat("  Optimizing grid point", pt, "\n")
  inp <- all_inputs[all_inputs$grid_point == pt, ]
  locest_results[idx, ] <- run_optimizer_on_arrays(
    inp$zi, inp$zj, inp$xl, inp$yl,
    df_true, alpha_true,
    ensemble = locest_ensemble, seed = 42, max_boundary_retries = 5,
    lower_bounds = c(0.01, 0.01), upper_bounds = c(upper_a, upper_b)
  )
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

save_csv(data.frame(
  grid_point = test_points,
  X = X[test_points],
  Y = Y[test_points],
  lower_a = rep(0.01, length(test_points)),
  lower_b = rep(0.01, length(test_points)),
  lower_g = rep(-pi / 2, length(test_points)),
  upper_a = rep(upper_a, length(test_points)),
  upper_b = rep(upper_b, length(test_points)),
  upper_g = rep(pi / 2, length(test_points)),
  ensemble = rep(locest_ensemble, length(test_points)),
  abstand = rep(locest_abst, length(test_points)),
  a_est = locest_results[, 1],
  b_est = locest_results[, 2],
  g_est = locest_results[, 3]
), "local_optimizer_selected_outputs",
"Exact R-selected optimizer outputs and bounds for the 5 selected local-estimation cells")

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
# Read back from CSV so downstream hclust uses exactly the serialized values
ell_dist <- as.matrix(read.csv(file.path(outdir, "ellipse_dissimilarity_matrix.csv")))

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

# ── 14c. In-cluster re-estimation (LEC) ──────────────────────────────────
cat("\n[Step 14c] In-cluster Re-estimation (LEC)\n")

# Read back simulation data from binary to avoid R vs Python CSV parser differences
sim_bin_path <- file.path(outdir, "simulation_nonstat_seed42.bin")
con <- file(sim_bin_path, "rb")
sim_bin_flat <- matrix(NA, nrow = n_grid, ncol = n_sim)
for (k in 1:n_sim) {
  sim_bin_flat[, k] <- readBin(con, what = "double", n = n_grid, size = 8, endian = "little")
}
close(con)
# Reconstruct 3D array from binary (column-major flat order, matching Python test)
sim_3d_rb <- array(NA, dim = c(nrow_grid, ncol_grid, n_sim))
for (r_idx in 1:n_grid) {
  r_i <- ((r_idx - 1) %% nrow_grid) + 1
  r_j <- ((r_idx - 1) %/% nrow_grid) + 1
  for (k in 1:n_sim) {
    sim_3d_rb[r_i, r_j, k] <- sim_bin_flat[r_idx, k]
  }
}

upper_bounds_locest <- c(max(a_matrix) + 5, 2 * max(b_matrix))
max_cl <- max(clusters_lec)
incluster_lec <- matrix(-Inf, nrow = max_cl, ncol = 5)
max_dist_val <- 4 * ((max(X) - min(X)) / (resolution - 1))

for (cl in 1:max_cl) {
  which_cl <- which(clusters_lec == cl)
  if (length(which_cl) < 5) next

  cat("  Cluster", cl, ":", length(which_cl), "cells\n")

  # Build z matrix from readback sim data (matching Python's extraction)
  z <- t(sapply(which_cl, function(j) {
    sim_3d_rb[number_grid(j)[1], number_grid(j)[2], ]
  }))

  n_cl <- length(which_cl)
  X_cl <- X[which_cl]
  Y_cl <- Y[which_cl]

  # Build pair arrays matching Python's pairwise_density_optim
  # Python: np.triu_indices(n_grid, k=1) → same order as R's lower.tri
  ilist <- rep(1, n_cl) %*% t(1:n_cl)
  ilist <- ilist[lower.tri(ilist)]
  jlist <- (1:n_cl) %*% t(rep(1, n_cl))
  jlist <- jlist[lower.tri(jlist)]

  dx <- X_cl[ilist] - X_cl[jlist]
  dy <- Y_cl[ilist] - Y_cl[jlist]

  # Filter by max_dist
  sel <- which(dx * dx + dy * dy <= max_dist_val * max_dist_val)
  if (length(sel) == 0) next

  ilist <- ilist[sel]
  jlist <- jlist[sel]
  dx <- dx[sel]
  dy <- dy[sel]

  # Expand by n_sim (matching Python's np.repeat and z[ilist].reshape(-1))
  zilist <- as.vector(t(z[ilist, ]))
  zjlist <- as.vector(t(z[jlist, ]))
  Xlist <- rep(dx, each = n_sim)
  Ylist <- rep(dy, each = n_sim)

  # Run optimizer matching Python's _optimize_local_mle_via_r with max_boundary_retries=0
  est <- run_optimizer_on_arrays(
    zilist, zjlist, Xlist, Ylist,
    df_true, alpha_true,
    ensemble = 3, seed = 42, max_boundary_retries = 0,
    lower_bounds = c(0.01, 0.01), upper_bounds = upper_bounds_locest
  )

  incluster_lec[cl, 1:3] <- est
  incluster_lec[cl, 4] <- length(which_cl)

  # Compute llh_in_cluster from the same readback pair arrays
  lh <- sum(pairwise_density_summand(zilist, zjlist, Xlist, Ylist,
                                     df_true, alpha_true, est[1], est[2], est[3]))
  incluster_lec[cl, 5] <- lh / length(zilist)
}

incluster_lec_df <- data.frame(
  cluster = 1:nrow(incluster_lec),
  a = incluster_lec[, 1],
  b = incluster_lec[, 2],
  g = incluster_lec[, 3],
  n_cells = incluster_lec[, 4],
  avg_llh = incluster_lec[, 5]
)
save_csv(incluster_lec_df, "calc_estimates_in_clusters_lec",
         "In-cluster re-estimation results for LEC clusters")

# ── 14b. crop_matrix / crop_local_estimates ───────────────────────────────
cat("\n[Step 14b] crop_matrix / crop_local_estimates\n")

for (margin in c(0, 1, 2)) {
  cropped <- crop_local_estimates(smoothed, margin)
  save_csv(as.data.frame(cropped),
           paste0("crop_local_estimates_margin", margin),
           paste0("crop_local_estimates(smoothed, margin=", margin, ")"))
}

# Also save crop_matrix on the a-column for direct testing
for (margin in c(1, 2)) {
  cropped_a <- as.vector(crop_matrix(smoothed[,1], margin))
  save_csv(data.frame(a_cropped = cropped_a),
           paste0("crop_matrix_a_margin", margin),
           paste0("crop_matrix(smoothed[,1], margin=", margin, ")"))
}

# ── 15.  Summary ──────────────────────────────────────────────────────────
cat("\n[Step 15] Summary\n")

# dist_x / dist_y helper test cases
dist_tests <- data.frame(
  x1 = c(-5.0, 0.0, 3.5, -1.25, 8.0),
  x2 = c(2.0, 0.0, 1.25, 4.75, -3.0),
  y1 = c(5.0, -3.0, 0.5, 2.25, -7.5),
  y2 = c(-1.0, 1.25, 0.5, -4.75, 2.5)
)
dist_tests$dist_x <- mapply(dist_x, dist_tests$x1, dist_tests$x2)
dist_tests$dist_y <- mapply(dist_y, dist_tests$y1, dist_tests$y2)
save_csv(dist_tests, "dist_helper_test_cases",
         "dist_x(x1,x2) and dist_y(y1,y2) helper results")

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

library(lhs)

source("r_code/functions.R")
source("r_code/parameters.R")

sim_ref <- read.csv("tests/reference_data/simulation_nonstat_seed42.csv")
sim_cols <- grep("^sim", names(sim_ref), value = TRUE)

a_sim <- array(NA_real_, dim = c(nrow_grid, ncol_grid, length(sim_cols)))
for (idx in 1:n_grid) {
  ij <- number_grid(idx)
  for (k in seq_along(sim_cols)) {
    a_sim[ij[1], ij[2], k] <- sim_ref[idx, sim_cols[k]]
  }
}

inspect_point <- function(pt) {
  x <- X[pt]
  y <- Y[pt]
  xy_pos <- which.min(dist_x(x, X)^2 + dist_y(y, Y)^2)[1]

  xs <- rep(-locest_abst:locest_abst, (2 * locest_abst + 1))
  ys <- rep(-locest_abst:locest_abst, each = (2 * locest_abst + 1))
  which.dist <- which(sqrt(xs * xs + ys * ys) <= locest_abst & (xs * xs + ys * ys) > 0)
  xs <- xs[which.dist]
  ys <- ys[which.dist]

  sel.x <- number_grid(xy_pos)[2] + xs
  sel.y <- number_grid(xy_pos)[1] + ys
  in.bounds <- which(sel.y > 0 & sel.y <= nrow_grid & sel.x > 0 & sel.x <= ncol_grid)
  sel.x <- sel.x[in.bounds]
  sel.y <- sel.y[in.bounds]
  sel.grid <- sapply(1:length(in.bounds), function(i) grid_number(sel.y[i], sel.x[i]))

  n_sel <- length(sel.x)
  zlist1 <- as.vector(sapply(1:n_sel, function(i) {
    a_sim[number_grid(sel.grid[i])[1], number_grid(sel.grid[i])[2], ]
  }))
  zlist2 <- rep(a_sim[number_grid(xy_pos)[1], number_grid(xy_pos)[2], ], n_sel)
  Xlist <- rep(X[sel.grid] - X[xy_pos], each = n_sim)
  Ylist <- rep(Y[sel.grid] - Y[xy_pos], each = n_sim)

  llh <- function(par) {
    sum(pairwise_density_summand(zlist1, zlist2, Xlist, Ylist, df_true, alpha_true,
                                 par[1], par[2], par[3]))
  }

  lo <- c(0.01, 0.01, -pi / 2)
  hi <- c(max(a_matrix) + 5, 2 * max(b_matrix), pi / 2)
  parscale <- (hi - lo) / 100

  set.seed(42)
  starts <- matrix(lo, locest_ensemble + 5, 3, byrow = TRUE) +
    maximinLHS(locest_ensemble + 5, 3) %*% diag(hi - lo)

  cat("\nPOINT", pt, "x=", x, "y=", y, "\n")
  print(round(starts, 6))

  num_calc <- 1
  num_more <- 0
  best <- NULL
  while (num_calc <= locest_ensemble) {
    start <- starts[num_calc + num_more, ]
    fit <- optim(start, fn = llh, method = "L-BFGS-B",
                 lower = lo, upper = hi,
                 control = list(fnscale = -1, parscale = parscale, maxit = 10000))
    wrapped <- FALSE
    if (abs(fit$par[3]) == pi / 2) {
      fit <- optim(c(fit$par[1], fit$par[2], -fit$par[3]), fn = llh, method = "L-BFGS-B",
                   lower = lo, upper = hi,
                   control = list(fnscale = -1, parscale = parscale, maxit = 10000))
      wrapped <- TRUE
    }
    near_bound <- min(abs(fit$par - lo)) < 0.01 || min(abs(fit$par - hi)) < 0.01
    cat(sprintf(
      "run=%d extra=%d start=(%.6f, %.6f, %.6f) fit=(%.6f, %.6f, %.6f) val=%.10f wrapped=%s near_bound=%s\n",
      num_calc, num_more, start[1], start[2], start[3],
      fit$par[1], fit$par[2], fit$par[3], fit$value, wrapped, near_bound
    ))
    if (is.null(best) || (-fit$value < -best$value)) {
      best <- fit
    }
    if (near_bound && num_more < 5) {
      num_more <- num_more + 1
      next
    }
    num_calc <- num_calc + 1
  }

  cat(sprintf(
    "BEST point=%d par=(%.12f, %.12f, %.12f) value=%.12f\n",
    pt, best$par[1], best$par[2], best$par[3], best$value
  ))
}

inspect_point(1)
inspect_point(100)
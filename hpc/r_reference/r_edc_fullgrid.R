#!/usr/bin/env Rscript
###############################################################################
# R reference for criterion F1: the deterministic EDC chain on one input.
#
# Chain (the same steps as tests/generate_cmip6_mini_reference.R):
#   monthly pr -> stl(s.window = "periodic", robust = TRUE) trend removed
#   -> annual maxima of the detrended months -> GEV fit (optim Nelder-Mead)
#   -> unit Frechet -> c_extrcoeff_matrix(madogram = TRUE) from
#   r_code/functions.R -> empirical quantile of the upper triangle
#   -> hclust(average) via clustering() -> cluster_number_threshold_method
#   -> max(k, min_k) -> cutree.
#
# The madogram, threshold, k rule and hclust call are Thorsten Dickhaus's
# functions, sourced unchanged from r_code/functions.R.  STL, annual maxima,
# GEV and Frechet steps are copied from the frozen mini-reference generator.
#
# Input modes (exactly one):
#   --netcdf-dir DIR    CMIP6 Amon pr files, read with ncdump (no R packages).
#                       Months are selected by decoded time in
#                       [--year-start, --year-end].
#   --monthly-csv FILE  (n_months x n_cells) CSV with header, cells in R
#                       column-major order, months consecutive from January.
#                       Needs --nrow and --ncol.
#   --monthly-bin FILE  the same matrix as little-endian float64, months
#                       fastest.  Needs --n-months, --nrow and --ncol.
#   --sim-csv FILE      CSV with flat_index and one column per replicate,
#                       cells in R column-major order.  Skips STL, maxima,
#                       GEV and Frechet; the replicates go straight into
#                       c_extrcoeff_matrix.  Needs --nrow and --ncol.
#   --sim-bin FILE      the same as little-endian float64, n_cells values per
#                       replicate.  Needs --n-sim, --nrow and --ncol.
#   --madogram-csv FILE a full madogram matrix read with read.csv, then only
#                       the threshold and clustering steps.  This replays the
#                       frozen generator, which read its matrix back this way.
#   --dump-input FILE   reader check: after loading, write the monthly input
#                       (months x lat x lon, f64le, months fastest) and stop.
#
# Exactness of text input.  R's decimal parser is not correctly rounded on
# every platform.  On aarch64 macOS (long double = double) only 83% of
# %.17g strings parse back to the double that was printed; on albedo's
# x86_64 R all did in a 1e5-value test.  Use the -bin modes when the input
# must be bit-exact.  Outputs are printed with %.17e and are exact for any
# correctly rounded parser (Python, NumPy).
#
# Cell orders.  "r" is R's column-major order over [lat, lon]
# (r_index = (lon_idx - 1) * n_lat + lat_idx).  "py" is weatherisk's
# row-major order (py_index = (lat_idx - 1) * n_lon + (lon_idx - 1), 0-based).
# hclust breaks ties by position, so both orders are clustered and compared.
#
# Outputs in --out (all floats written with %.17e, 18 significant digits):
#   run_info.csv                 parameters, sizes, timings, diagnostics
#   grid_cells.csv               r_index, py_index, lat_idx, lon_idx, lat, lon
#   gev_params.csv               per cell loc, scale, shape, convergence, nll
#   stl_trend.f64le.bin.gz       trend, n_months x n_cells, R cell order,
#                                months fastest (column-major), little endian
#   annual_maxima.csv.gz         one row per cell, one column per year
#   frechet.csv.gz               one row per cell, one column per year
#   madogram_sum_condensed_py.u16le.bin.gz
#                                S = sum_k |r_ik - r_jk| per pair, uint16,
#                                scipy condensed order in py cell order
#   madogram_value_lut.csv       S -> the exact R double madogram value
#   hclust_<order>_merge.csv     merge_i, merge_j, height (R hclust format)
#   hclust_<order>_order.csv     leaf order
#   edc_summary.csv              per order q, k_raw, k, ties at threshold
#   labels.csv                   per cell labels in both orders
#   madogram_matrix.csv          only with --save-matrix (small grids)
#   stl_trend.csv, detrended.csv only with --save-matrix (small grids)
###############################################################################

t_start <- proc.time()[["elapsed"]]

# ── Arguments ──────────────────────────────────────────────────────────────
parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    key <- sub("^--", "", a)
    if (grepl("=", key, fixed = TRUE)) {
      out[[sub("=.*$", "", key)]] <- sub("^[^=]*=", "", key)
      i <- i + 1
    } else if (i < length(argv) && !startsWith(argv[i + 1], "--")) {
      out[[key]] <- argv[i + 1]
      i <- i + 2
    } else {
      out[[key]] <- "TRUE"
      i <- i + 1
    }
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(name, default = NULL) {
  if (!is.null(opt[[name]])) opt[[name]] else default
}

functions_path <- get_opt("functions", "r_code/functions.R")
out_dir        <- get_opt("out")
quantile_level <- as.numeric(get_opt("quantile", "0.3"))
min_k          <- as.integer(get_opt("min-k", "2"))
period         <- as.integer(get_opt("period", "12"))
year_start     <- as.integer(get_opt("year-start", "1850"))
year_end       <- as.integer(get_opt("year-end", "2005"))
ncdump_bin     <- get_opt("ncdump", "ncdump")
orders         <- strsplit(get_opt("orders", "r,py"), ",", fixed = TRUE)[[1]]
save_matrix    <- identical(get_opt("save-matrix", "FALSE"), "TRUE")

modes <- c("netcdf-dir", "monthly-csv", "monthly-bin", "sim-csv", "sim-bin", "madogram-csv")
given <- modes[vapply(modes, function(m) !is.null(opt[[m]]), logical(1))]
if (length(given) != 1) stop("Give exactly one of --", paste(modes, collapse = ", --"))
if (is.null(out_dir)) stop("--out is required")
if (!all(orders %in% c("r", "py"))) stop("--orders must be a subset of r,py")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

info <- list()
add_info <- function(key, value) {
  info[[key]] <<- if (is.numeric(value)) sprintf("%.17g", value) else as.character(value)
}
fmt <- function(x) sprintf("%.17e", x)
logmsg <- function(...) {
  cat(sprintf("[%7.1fs] ", proc.time()[["elapsed"]] - t_start), ..., "\n", sep = "")
  flush.console()
}

source(functions_path)   # Thorsten's functions; defines no globals we rely on

# ── Input readers ──────────────────────────────────────────────────────────
run_capture <- function(cmd, args) {
  res <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    stop(cmd, " ", paste(args, collapse = " "), " failed:\n", paste(res, collapse = "\n"))
  }
  res
}

# Values of one variable printed by ncdump, as character tokens.
ncdump_values <- function(file, var, digits) {
  cmd <- sprintf(
    "%s -v %s -p %s %s | awk '/^ %s =/{f=1; sub(/^ %s =/,\"\")} f{gsub(/[;}]/,\"\"); n=split($0,a,\",\"); for(i=1;i<=n;i++){g=a[i]; gsub(/[ \\t]/,\"\",g); if(g!=\"\") print g}}'",
    ncdump_bin, var, digits, shQuote(file), var, var)
  con <- pipe(cmd, open = "r")
  on.exit(close(con))
  scan(con, what = character(), quiet = TRUE)
}

# float32 values printed with 9 significant digits identify the float
# uniquely.  Rounding the parsed double to float32 recovers the stored value
# exactly, which is what xarray hands to weatherisk before the float64 cast.
to_float32 <- function(x) {
  readBin(writeBin(as.numeric(x), raw(), size = 4, endian = "little"),
          "double", n = length(x), size = 4, endian = "little")
}

read_netcdf_dir <- function(dir) {
  files <- sort(list.files(dir, pattern = "\\.nc$", full.names = TRUE))
  if (length(files) == 0) stop("No .nc files in ", dir)
  blocks <- list(); stamps <- character(); lat <- NULL; lon <- NULL
  n_files_used <- 0
  for (f in files) {
    hdr <- run_capture(ncdump_bin, c("-h", shQuote(f)))
    if (!any(grepl("^\\s*float pr\\(time, lat, lon\\)", hdr)))
      stop(f, ": pr is not float pr(time, lat, lon)")
    if (any(grepl("pr:(scale_factor|add_offset)", hdr)))
      stop(f, ": pr has scale_factor/add_offset; the reader does not apply them")
    tl <- run_capture(ncdump_bin, c("-t", "-v", "time", shQuote(f)))
    tl <- paste(tl[seq(which(grepl("^ time =", tl)), length(tl))], collapse = " ")
    stamp <- regmatches(tl, gregexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", tl))[[1]]
    yrs <- as.integer(substr(stamp, 1, 4))
    keep <- which(yrs >= year_start & yrs <= year_end)
    if (length(keep) == 0) next
    n_lat_f <- as.integer(sub(".*lat = ([0-9]+) ;.*", "\\1", grep("^\\s*lat = [0-9]+ ;", hdr, value = TRUE)))
    n_lon_f <- as.integer(sub(".*lon = ([0-9]+) ;.*", "\\1", grep("^\\s*lon = [0-9]+ ;", hdr, value = TRUE)))
    tok <- ncdump_values(f, "pr", "9")
    if (length(tok) != length(stamp) * n_lat_f * n_lon_f)
      stop(f, ": read ", length(tok), " values, expected ", length(stamp) * n_lat_f * n_lon_f)
    if (any(tok == "_")) stop(f, ": contains fill values; weatherisk would see NaN here")
    vals <- to_float32(as.numeric(tok))
    if (any(!is.finite(vals)) || any(abs(vals) >= 1e19)) stop(f, ": non-finite or fill-like values")
    dim(vals) <- c(n_lon_f, n_lat_f, length(stamp))      # C order: lon fastest
    blocks[[length(blocks) + 1]] <- aperm(vals, c(3, 2, 1))[keep, , , drop = FALSE]
    stamps <- c(stamps, stamp[keep])
    if (is.null(lat)) {
      lat <- as.numeric(ncdump_values(f, "lat", "17"))
      lon <- as.numeric(ncdump_values(f, "lon", "17"))
      if (length(lat) != n_lat_f || length(lon) != n_lon_f)
        stop(f, ": read ", length(lat), " lat and ", length(lon), " lon values, header says ",
             n_lat_f, " and ", n_lon_f)
    }
    n_files_used <- n_files_used + 1
    logmsg("read ", basename(f), " (", length(keep), " months)")
  }
  n_months <- sum(vapply(blocks, function(b) dim(b)[1], integer(1)))
  pr <- array(0, dim = c(n_months, length(lat), length(lon)))
  pos <- 0
  for (b in blocks) {
    pr[pos + seq_len(dim(b)[1]), , ] <- b
    pos <- pos + dim(b)[1]
  }
  ym <- as.integer(substr(stamps, 1, 4)) * 12 + as.integer(substr(stamps, 6, 7)) - 1
  if (any(diff(ym) != 1)) stop("Selected months are not consecutive")
  add_info("n_files_used", n_files_used)
  add_info("first_month", substr(stamps[1], 1, 7))
  add_info("last_month", substr(stamps[length(stamps)], 1, 7))
  list(pr = pr, year = as.integer(substr(stamps, 1, 4)), lat = lat, lon = lon)
}

read_monthly_csv <- function(path, nrow_g, ncol_g) {
  m <- as.matrix(read.csv(path))
  if (ncol(m) != nrow_g * ncol_g) stop("monthly CSV has ", ncol(m), " cells")
  if (nrow(m) %% period != 0) stop("monthly CSV rows not a multiple of the period")
  pr <- array(m, dim = c(nrow(m), nrow_g, ncol_g))
  list(pr = pr, year = rep(seq_len(nrow(m) / period), each = period),
       lat = seq_len(nrow_g), lon = seq_len(ncol_g))
}

# ── Load input ─────────────────────────────────────────────────────────────
mode <- given
add_info("mode", mode)
add_info("input", opt[[mode]])
add_info("functions", normalizePath(functions_path))
add_info("R_version", R.version.string)

read_f64 <- function(path, n) {
  con <- file(path, "rb"); on.exit(close(con))
  x <- readBin(con, "double", n = n + 1, size = 8, endian = "little")
  if (length(x) != n) stop(path, " holds ", length(x), " doubles, expected ", n)
  x
}
is_sim  <- mode %in% c("sim-csv", "sim-bin")
is_mado <- mode == "madogram-csv"
if (is_sim || is_mado) {
  nrow_g <- as.integer(get_opt("nrow")); ncol_g <- as.integer(get_opt("ncol"))
  lat_vals <- seq_len(nrow_g); lon_vals <- seq_len(ncol_g)
}
if (mode == "sim-csv") {
  sim_df <- read.csv(opt[["sim-csv"]])
  if (!is.null(sim_df$flat_index) && !identical(as.integer(sim_df$flat_index), seq_len(nrow(sim_df))))
    stop("sim CSV flat_index must be 1..n in R order")
  sim_mat <- as.matrix(sim_df[, setdiff(names(sim_df), "flat_index")])
} else if (mode == "sim-bin") {
  k_sim <- as.integer(get_opt("n-sim"))
  sim_mat <- matrix(read_f64(opt[["sim-bin"]], nrow_g * ncol_g * k_sim), nrow_g * ncol_g, k_sim)
} else if (is_mado) {
  v_matrix <- unname(as.matrix(read.csv(opt[["madogram-csv"]])))
  if (nrow(v_matrix) != nrow_g * ncol_g) stop("madogram CSV has ", nrow(v_matrix), " rows")
}
if (is_sim) {
  if (nrow(sim_mat) != nrow_g * ncol_g) stop("sim input has ", nrow(sim_mat), " cells")
  n_years <- ncol(sim_mat)
} else if (!is_mado) {
  if (mode == "netcdf-dir") {
    inp <- read_netcdf_dir(opt[["netcdf-dir"]])
  } else if (mode == "monthly-bin") {
    nm <- as.integer(get_opt("n-months"))
    nr <- as.integer(get_opt("nrow")); nc <- as.integer(get_opt("ncol"))
    inp <- list(pr = array(read_f64(opt[["monthly-bin"]], nm * nr * nc), dim = c(nm, nr, nc)),
                year = rep(seq_len(nm / period), each = period),
                lat = seq_len(nr), lon = seq_len(nc))
    if (nm %% period != 0) stop("--n-months must be a multiple of the period")
  } else {
    inp <- read_monthly_csv(opt[["monthly-csv"]], as.integer(get_opt("nrow")),
                            as.integer(get_opt("ncol")))
  }
  pr <- inp$pr; lat_vals <- inp$lat; lon_vals <- inp$lon
  nrow_g <- dim(pr)[2]; ncol_g <- dim(pr)[3]
  logmsg("input: ", dim(pr)[1], " months x ", nrow_g, " x ", ncol_g)
  if (!is.null(opt[["dump-input"]])) {
    # Reader check only: write the monthly input and stop.
    con <- file(opt[["dump-input"]], "wb")
    writeBin(as.vector(pr), con, size = 8, endian = "little"); close(con)
    write.csv(data.frame(parameter = names(info), value = unlist(info)),
              file.path(out_dir, "run_info.csv"), row.names = FALSE)
    logmsg("wrote ", opt[["dump-input"]], " (", paste(dim(pr), collapse = " x "), ", f64le, months fastest)")
    quit(status = 0)
  }
}

# Globals that r_code/functions.R reads.
nrow_grid <- nrow_g
ncol_grid <- ncol_g
n_grid    <- nrow_g * ncol_g
add_info("n_lat", nrow_g); add_info("n_lon", ncol_g); add_info("n_cells", n_grid)

r_lat <- rep(seq_len(nrow_g), ncol_g)            # lat_idx of R cell s
r_lon <- rep(seq_len(ncol_g), each = nrow_g)     # lon_idx of R cell s
py_of_r <- (r_lat - 1L) * ncol_g + (r_lon - 1L)  # 0-based weatherisk index
r_of_py <- order(py_of_r)                        # R cell at py position k+1

write.csv(data.frame(r_index = seq_len(n_grid), py_index = py_of_r,
                     lat_idx = r_lat, lon_idx = r_lon,
                     lat = fmt(lat_vals[r_lat]), lon = fmt(lon_vals[r_lon])),
          file.path(out_dir, "grid_cells.csv"), row.names = FALSE, quote = FALSE)

# ── STL, annual maxima, GEV, Frechet ───────────────────────────────────────
if (!is_sim && !is_mado) {
  n_months <- dim(pr)[1]
  trend <- matrix(0, n_months, n_grid)
  t0 <- proc.time()[["elapsed"]]
  for (s in seq_len(n_grid)) {
    x <- pr[, r_lat[s], r_lon[s]]
    fit <- stl(ts(x, frequency = period), s.window = "periodic", robust = TRUE)
    trend[, s] <- fit$time.series[, "trend"]
    if (s %% 2000 == 0) logmsg("stl ", s, "/", n_grid)
  }
  add_info("seconds_stl", proc.time()[["elapsed"]] - t0)
  pr_mat <- matrix(pr, n_months, n_grid)          # same R cell order
  detrended <- pr_mat - trend
  logmsg("stl done")

  con <- gzfile(file.path(out_dir, "stl_trend.f64le.bin.gz"), "wb")
  writeBin(as.vector(trend), con, size = 8, endian = "little")
  close(con)
  if (save_matrix) {
    write.csv(matrix(fmt(trend), n_months), file.path(out_dir, "stl_trend.csv"),
              row.names = FALSE, quote = FALSE)
    write.csv(matrix(fmt(detrended), n_months), file.path(out_dir, "detrended.csv"),
              row.names = FALSE, quote = FALSE)
  }

  years <- unique(inp$year)
  complete <- years[vapply(years, function(y) sum(inp$year == y) == period, logical(1))]
  if (length(complete) != length(years)) logmsg("dropping incomplete years: ",
                                                 paste(setdiff(years, complete), collapse = " "))
  n_years <- length(complete)
  am <- matrix(0, n_years, n_grid)
  for (k in seq_len(n_years)) {
    am[k, ] <- apply(detrended[inp$year == complete[k], , drop = FALSE], 2, max)
  }
  rm(pr, pr_mat, detrended); invisible(gc())
  add_info("n_years", n_years)
  add_info("year_first", complete[1]); add_info("year_last", complete[n_years])

  # weatherisk's validity rule: not all NaN and population sd >= 1e-15.
  pop_sd <- apply(am, 2, function(v) sqrt(mean((v - mean(v))^2)))
  n_invalid <- sum(!is.finite(pop_sd) | pop_sd < 1e-15)
  add_info("n_invalid_cells", n_invalid)
  if (n_invalid > 0) stop(n_invalid, " cells fail the validity rule; the R chain needs a full grid")

  gev <- data.frame(loc = NA_real_, scale = NA_real_, shape = NA_real_,
                    convergence = NA_integer_, nll = NA_real_)[rep(1, n_grid), ]
  frechet_mat <- am
  n_clipped <- 0
  t0 <- proc.time()[["elapsed"]]
  for (s in seq_len(n_grid)) {
    y_data <- am[, s]
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
    mu0    <- mean(y_data)
    sig0   <- sd(y_data) * sqrt(6) / pi
    init   <- c(mu0 - 0.5772 * sig0, sig0, 0.1)
    result <- optim(init, gev_nll, method = "Nelder-Mead",
                    control = list(maxit = 5000))
    mu <- result$par[1]; sig <- result$par[2]; xi <- result$par[3]
    gev[s, ] <- list(mu, sig, xi, result$convergence, result$value)
    if (abs(xi) < 1e-8) {
      frechet_mat[, s] <- exp((y_data - mu) / sig)
    } else {
      w <- 1 + xi * (y_data - mu) / sig
      n_clipped <- n_clipped + sum(w < 1e-10)
      w <- pmax(w, 1e-10)
      frechet_mat[, s] <- w^(1 / xi)
    }
  }
  add_info("seconds_gev", proc.time()[["elapsed"]] - t0)
  add_info("n_gev_not_converged", sum(gev$convergence != 0))
  logmsg("gev done")

  write.csv(data.frame(r_index = seq_len(n_grid), py_index = py_of_r,
                       loc = fmt(gev$loc), scale = fmt(gev$scale), shape = fmt(gev$shape),
                       convergence = gev$convergence, nll = fmt(gev$nll)),
            file.path(out_dir, "gev_params.csv"), row.names = FALSE, quote = FALSE)
  write_cells <- function(mat, name) {
    df <- data.frame(r_index = seq_len(n_grid), py_index = py_of_r,
                     matrix(fmt(t(mat)), n_grid))
    names(df)[-(1:2)] <- paste0("y", complete)
    con <- gzfile(file.path(out_dir, name), "w")
    write.csv(df, con, row.names = FALSE, quote = FALSE)
    close(con)
  }
  write_cells(am, "annual_maxima.csv.gz")
  write_cells(frechet_mat, "frechet.csv.gz")

  # Tie diagnostics.  The madogram depends only on within-cell ranks, so a
  # tie changes it only through the tie convention (R counts values >= x,
  # weatherisk averages).  Clipping can create ties absent from the maxima.
  r_rank <- function(v) vapply(v, function(a) sum(a <= v), numeric(1))
  n_tie_am <- sum(apply(am, 2, function(v) anyDuplicated(v) > 0))
  n_tie_fr <- sum(apply(frechet_mat, 2, function(v) anyDuplicated(v) > 0))
  n_rank_changed <- sum(vapply(seq_len(n_grid), function(s)
    !identical(r_rank(am[, s]), r_rank(frechet_mat[, s])), logical(1)))
  add_info("n_cells_with_ties_in_annual_maxima", n_tie_am)
  add_info("n_cells_with_ties_in_frechet", n_tie_fr)
  add_info("n_cells_rank_changed_by_frechet", n_rank_changed)
  add_info("n_cells_frechet_below_0.05", sum(apply(frechet_mat, 2, min) < 0.05))
  add_info("n_cells_frechet_above_nyears_sq", sum(apply(frechet_mat, 2, max) > n_years^2))
  add_info("n_frechet_values_clipped_by_r_floor", n_clipped)
  rm(am); invisible(gc())
} else if (is_sim) {
  frechet_mat <- t(sim_mat)                       # n_sim x n_grid
  n_years <- ncol(sim_mat)
  add_info("n_years", n_years)
  add_info("n_cells_with_ties_in_frechet",
           sum(apply(frechet_mat, 2, function(v) anyDuplicated(v) > 0)))
}

# ── Madogram (Thorsten's c_extrcoeff_matrix, unchanged) ────────────────────
if (!is_mado) {
  n_sim <- n_years
  a_sim <- array(t(frechet_mat), dim = c(nrow_g, ncol_g, n_years))
  rm(frechet_mat); invisible(gc())
  t0 <- proc.time()[["elapsed"]]
  v_matrix <- c_extrcoeff_matrix(a_sim, madogram = TRUE)
  add_info("seconds_madogram", proc.time()[["elapsed"]] - t0)
  rm(a_sim); invisible(gc())
  logmsg("madogram done")
}

if (save_matrix) {
  write.csv(matrix(fmt(v_matrix), n_grid), file.path(out_dir, "madogram_matrix.csv"),
            row.names = FALSE, quote = FALSE)
}

v_upper <- v_matrix[upper.tri(v_matrix)]
q_thr <- unname(quantile(v_upper, quantile_level))
rm(v_upper); invisible(gc())
add_info("quantile_level", quantile_level)
add_info("q_threshold", q_thr)
logmsg("threshold ", fmt(q_thr))

# ── Clustering in each cell order ──────────────────────────────────────────
summary_rows <- list()
labels <- data.frame(r_index = seq_len(n_grid), py_index = py_of_r,
                     lat_idx = r_lat, lon_idx = r_lon)
heights <- list()
for (ord in orders) {
  t0 <- proc.time()[["elapsed"]]
  if (ord == "r") {
    hc <- clustering(v_matrix)                     # Thorsten's call
  } else {
    d_py <- as.dist(v_matrix[r_of_py, r_of_py])    # scipy condensed order
  }
  if (ord == "py" && !is_mado) {
    # Exact integer numerator S of every madogram value, plus the map
    # S -> R double.  2 bytes per pair instead of 8, and lossless.
    denom <- n_sim * 2 * (n_sim + 1)
    S <- as.integer(round(as.vector(d_py) * denom))
    lut <- rep(NA_real_, max(S) + 1)
    lut[S + 1] <- as.vector(d_py)
    if (!all(lut[S + 1] == as.vector(d_py))) stop("S -> madogram value is not one-to-one")
    con <- gzfile(file.path(out_dir, "madogram_sum_condensed_py.u16le.bin.gz"), "wb")
    chunk <- 5e7
    for (st in seq(1, length(S), by = chunk)) {
      writeBin(S[st:min(length(S), st + chunk - 1)], con, size = 2, endian = "little")
    }
    close(con)
    used <- which(!is.na(lut))
    write.csv(data.frame(S = used - 1, value = fmt(lut[used])),
              file.path(out_dir, "madogram_value_lut.csv"), row.names = FALSE, quote = FALSE)
    add_info("madogram_denominator", denom)
    add_info("n_pairs", length(S))
    rm(S); invisible(gc())
  }
  if (ord == "py") {
    # clustering() is hclust(as.dist(t(m), diag = TRUE), "average"); for a
    # symmetric m this is hclust(as.dist(m), "average").
    hc <- hclust(d_py, method = "average")
    rm(d_py); invisible(gc())
  }
  k_raw <- cluster_number_threshold_method(hc, q_thr)
  k <- max(k_raw, min_k)
  lab <- cutree(hc, k)
  if (ord == "py") lab <- lab[order(r_of_py)]       # back to R cell order
  labels[[paste0("label_", ord)]] <- as.integer(lab)
  heights[[ord]] <- hc$height
  write.csv(data.frame(merge_i = hc$merge[, 1], merge_j = hc$merge[, 2],
                       height = fmt(hc$height)),
            file.path(out_dir, paste0("hclust_", ord, "_merge.csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(data.frame(order = hc$order),
            file.path(out_dir, paste0("hclust_", ord, "_order.csv")),
            row.names = FALSE, quote = FALSE)
  summary_rows[[ord]] <- data.frame(
    order = ord, q_threshold = fmt(q_thr), k_raw = k_raw, k = k, min_k = min_k,
    n_heights_equal_threshold = sum(hc$height == q_thr),
    n_heights_within_1e12 = sum(abs(hc$height - q_thr) <= 1e-12),
    n_unique_heights = length(unique(hc$height)),
    seconds_hclust = round(proc.time()[["elapsed"]] - t0, 3))
  logmsg("order ", ord, ": k_raw=", k_raw, " k=", k)
  rm(hc); invisible(gc())
}
rm(v_matrix); invisible(gc())

summ <- do.call(rbind, summary_rows)
if (all(c("r", "py") %in% orders)) {
  ari <- function(x, y) {
    tab <- table(x, y); c2 <- function(n) n * (n - 1) / 2
    a <- sum(c2(tab)); b <- sum(c2(rowSums(tab))); c <- sum(c2(colSums(tab)))
    e <- b * c / c2(length(x)); (a - e) / ((b + c) / 2 - e)
  }
  add_info("ari_r_vs_py", ari(labels$label_r, labels$label_py))
  add_info("max_abs_sorted_height_diff_r_vs_py",
           max(abs(sort(heights$r) - sort(heights$py))))
}
write.csv(summ, file.path(out_dir, "edc_summary.csv"), row.names = FALSE, quote = FALSE)
write.csv(labels, file.path(out_dir, "labels.csv"), row.names = FALSE, quote = FALSE)

add_info("seconds_total", proc.time()[["elapsed"]] - t_start)
write.csv(data.frame(parameter = names(info), value = unlist(info)),
          file.path(out_dir, "run_info.csv"), row.names = FALSE)
print(summ)
logmsg("done")

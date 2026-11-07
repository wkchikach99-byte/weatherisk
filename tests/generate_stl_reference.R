#!/usr/bin/env Rscript
###############################################################################
# Generate frozen R references for weatherisk.stl (tests/test_stl_r_parity.py).
#
# Each case is a batch of synthetic series run through stats::stl column by
# column. Values are stored as raw little-endian float64, column-major, so no
# decimal rounding enters the comparison:
#
#   <case>.f64   x, seasonal, trend, remainder, weights   (5 * n * m doubles)
#   cases.csv    the stl() arguments plus R's effective win/deg/jump/inner/outer
#   meta.csv     R version and platform that produced the files
#
# Usage:  cd weatherisk && Rscript tests/generate_stl_reference.R
###############################################################################

outdir <- "tests/reference_data/stl"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261107)

gen <- function(n, m, period = 12, outliers = TRUE, dry = FALSE) {
  t <- 1:n
  sapply(1:m, function(k) {
    x <- 15 + runif(1, 1, 5) * sin(2 * pi * t / period + runif(1, 0, 6.3)) +
      runif(1, -3e-3, 3e-3) * t + rnorm(n, 0, 1.5)
    # Exact zeros make ties for psort and windows whose robustness weights
    # are all zero, which exercises the ok = FALSE branch of est.
    if (dry) x <- pmax(x - 16, 0)
    if (outliers) {
      o <- runif(n) < 0.03
      x <- x + o * sample(c(-1, 1), n, TRUE) * runif(n, 10, 20)
    }
    x
  })
}

cases <- list(
  list(name = "robust240",    x = gen(240, 25)),
  list(name = "robust1872",   x = gen(1872, 4)),
  list(name = "odd239",       x = gen(239, 3)),
  list(name = "odd241",       x = gen(241, 3)),
  list(name = "even250",      x = gen(250, 3)),
  list(name = "odd1871",      x = gen(1871, 2)),
  list(name = "dry240",       x = gen(240, 15, outliers = FALSE, dry = TRUE)),
  list(name = "nonrobust240", x = gen(240, 5), robust = FALSE),
  list(name = "swindow7",     x = gen(240, 3), s.window = 7),
  list(name = "sdegree1",     x = gen(240, 3), s.window = 13, s.degree = 1),
  list(name = "jump1",        x = gen(240, 3), jump1 = TRUE),
  list(name = "period4",      x = gen(120, 3, period = 4), period = 4),
  list(name = "period7",      x = gen(211, 3, period = 7), period = 7),
  list(name = "inner2outer5", x = gen(240, 3), inner = 2, outer = 5)
)

rows <- list()
for (cs in cases) {
  period <- if (is.null(cs$period)) 12 else cs$period
  args <- list(s.window = if (is.null(cs$s.window)) "periodic" else cs$s.window,
               robust = if (is.null(cs$robust)) TRUE else cs$robust)
  if (!is.null(cs$s.degree)) args$s.degree <- cs$s.degree
  if (isTRUE(cs$jump1)) { args$s.jump <- 1; args$t.jump <- 1; args$l.jump <- 1 }
  if (!is.null(cs$inner)) args$inner <- cs$inner
  if (!is.null(cs$outer)) args$outer <- cs$outer
  n <- nrow(cs$x); m <- ncol(cs$x)
  out <- array(0, c(n, m, 4))
  for (k in 1:m) {
    f <- do.call(stl, c(list(ts(cs$x[, k], frequency = period)), args))
    out[, k, ] <- cbind(f$time.series[, c("seasonal", "trend", "remainder")],
                        f$weights)
  }
  writeBin(c(as.vector(cs$x), as.vector(out)),
           file.path(outdir, paste0(cs$name, ".f64")), endian = "little")
  rows[[cs$name]] <- data.frame(
    case = cs$name, n = n, m = m, period = period,
    s_window = as.character(args$s.window),
    s_degree_arg = if (is.null(cs$s.degree)) NA else cs$s.degree,
    robust = args$robust, jump1 = isTRUE(cs$jump1),
    inner_arg = if (is.null(cs$inner)) NA else cs$inner,
    outer_arg = if (is.null(cs$outer)) NA else cs$outer,
    win_s = f$win[["s"]], win_t = f$win[["t"]], win_l = f$win[["l"]],
    deg_s = f$deg[["s"]], deg_t = f$deg[["t"]], deg_l = f$deg[["l"]],
    jump_s = f$jump[["s"]], jump_t = f$jump[["t"]], jump_l = f$jump[["l"]],
    inner = f$inner, outer = f$outer)
}
write.csv(do.call(rbind, rows), file.path(outdir, "cases.csv"), row.names = FALSE)
write.csv(data.frame(r_version = R.version.string, platform = R.version$platform),
          file.path(outdir, "meta.csv"), row.names = FALSE)
cat("Wrote", length(cases), "cases to", outdir, "\n")

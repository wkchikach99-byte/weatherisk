#!/usr/bin/env Rscript
###############################################################################
# Validate r_edc_fullgrid.R and r_lec_on_estimates.R against the frozen R
# references in tests/reference_data/ (read only).
#
# Runs each script through Rscript on frozen inputs, into a scratch
# directory, and compares the written outputs with the frozen files.
#
# Exactness.  The top-level references are printed with %.17e and the scripts
# print with %.17e, so two printed values are the same double if and only if
# the strings are equal.  Those checks compare strings and never parse.
# R's decimal parser is not correctly rounded on aarch64 macOS, so exact
# inputs are fed as binary.  One frozen CSV input (local_estimates_smoothed)
# has no binary twin; it is converted with system python3 (float() is
# correctly rounded; standard library only).
#
# The cmip6_mini references carry 15 significant digits (write.csv default),
# so their continuous stages are compared with a tolerance.
#
# Usage (from the repository root):
#   Rscript hpc/r_reference/validate_against_frozen.R SCRATCH_DIR
# Exit status 0 only if every check passes.
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: validate_against_frozen.R SCRATCH_DIR")
scratch <- args[1]
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
ref <- "tests/reference_data"
mini <- file.path(ref, "cmip6_mini")
here <- "hpc/r_reference"
fn <- c("--functions", "r_code/functions.R")

results <- list()
check <- function(name, ok, detail = "") {
  results[[length(results) + 1]] <<- data.frame(check = name, pass = isTRUE(ok), detail = detail)
  cat(sprintf("%-4s %-58s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}
note <- function(...) cat("     ", ..., "\n", sep = "")
num_cmp <- function(name, got, want, rtol) {
  got <- as.numeric(as.matrix(got)); want <- as.numeric(as.matrix(want))
  if (length(got) != length(want)) return(check(name, FALSE, sprintf("length %d vs %d", length(got), length(want))))
  tol <- rtol * max(1, max(abs(want)))
  d <- max(abs(got - want))
  check(name, d <= tol, sprintf("max abs diff %.3g (tol %.1e, n=%d)", d, tol, length(got)))
}
exact_cmp <- function(name, got, want) {
  got <- as.vector(as.matrix(got)); want <- as.vector(as.matrix(want))
  n_diff <- if (length(got) == length(want)) sum(got != want) else NA
  check(name, isTRUE(n_diff == 0), sprintf("%d values, %s differ", length(want), n_diff))
}
# Character tokens of a CSV with numeric fields (header dropped).  Integer
# tokens such as "0" or "100" are rewritten in %.17e so that they compare
# against %.17e output; an integer string parses exactly.
tokens <- function(path, cols = NULL) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  lines <- readLines(con); close(con)
  hdr <- gsub('"', "", strsplit(lines[1], ",", fixed = TRUE)[[1]])
  m <- do.call(rbind, strsplit(lines[-1], ",", fixed = TRUE))
  colnames(m) <- hdr
  if (!is.null(cols)) m <- m[, cols, drop = FALSE]
  int <- grepl("^-?[0-9]+$", m)
  m[int] <- sprintf("%.17e", as.numeric(m[int]))
  m
}
tok_cmp <- function(name, got, want) {
  got <- as.vector(got); want <- as.vector(want)
  n_diff <- if (length(got) == length(want)) sum(got != want) else NA
  check(name, isTRUE(n_diff == 0), sprintf("%d values, %s differ (string compare)", length(want), n_diff))
}
run <- function(script, argv) {
  log <- file.path(scratch, paste0(basename(script), ".log"))
  status <- system2("Rscript", c(file.path(here, script), argv), stdout = log, stderr = log)
  if (status != 0) stop(script, " failed (status ", status, "); see ", log)
}
rd <- function(...) read.csv(file.path(...))
cells_by_year <- function(path) {          # r_index, py_index, y... -> years x cells
  df <- read.csv(gzfile(path))
  t(as.matrix(df[order(df$r_index), -(1:2)]))
}
ari <- function(x, y) {
  tab <- table(x, y); c2 <- function(n) n * (n - 1) / 2
  a <- sum(c2(tab)); b <- sum(c2(rowSums(tab))); c <- sum(c2(colSums(tab)))
  e <- b * c / c2(length(x)); (a - e) / ((b + c) / 2 - e)
}
rel15 <- 5e-15     # half a unit in the 15th significant digit, plus slack
# The frozen references were generated with R on aarch64 macOS.  On another
# platform stl (Fortran), optim and libm round differently.  Measured on
# albedo (x86_64, R 4.2.2): STL trend 5.9e-13, Frechet 1.9e-11 absolute,
# i.e. at most 1.2e-13 relative to the largest value.  Continuous stages
# then use 1e-12 relative; every discrete result stays exact.
gen_platform <- R.version$arch == "aarch64" && Sys.info()[["sysname"]] == "Darwin"
rel_cont <- if (gen_platform) rel15 else 1e-12
cat(sprintf("platform %s %s: continuous-stage tolerance %.0e relative\n",
            Sys.info()[["sysname"]], R.version$arch, rel_cont))
# Does read.csv return exactly the doubles a %.17e file prints?  The frozen
# hclust references were computed on read.csv'd matrices, so they can be
# replayed only by a parser that rounds the same way as the generator's.
parser_exact <- local({
  t <- tokens(file.path(ref, "madogram_matrix.csv"))
  all(sprintf("%.17e", as.matrix(read.csv(file.path(ref, "madogram_matrix.csv")))) == as.vector(t))
})
cat(sprintf("read.csv round-trips %%.17e exactly: %s\n", parser_exact))

# ── 1. Mini CMIP6 EDC chain ────────────────────────────────────────────────
cat("\n== 1. Mini CMIP6 EDC chain (tests/reference_data/cmip6_mini) ==\n")
# The frozen monthly_pr.csv is rounded to 15 digits.  The generator's input
# is rebuilt exactly with its own code (step 1 of
# tests/generate_cmip6_mini_reference.R, copied unchanged) and passed as
# binary, so the chain runs on the generator's exact input.
resolution <- 5; n_years <- 20; nrow_grid <- 5; ncol_grid <- 5
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
num_cmp("mini rebuilt input == frozen monthly_pr.csv", matrix(pr, n_months), rd(mini, "monthly_pr.csv"), rel15)
pr_bin <- file.path(scratch, "mini_monthly_pr.f64le.bin")
con <- file(pr_bin, "wb"); writeBin(as.vector(pr), con, size = 8, endian = "little"); close(con)

o <- file.path(scratch, "mini_edc")
run("r_edc_fullgrid.R", c(fn, "--monthly-bin", pr_bin, "--n-months", "240", "--nrow", "5", "--ncol", "5",
                          "--out", o, "--save-matrix", "--min-k", "2"))
num_cmp("mini stl trend", rd(o, "stl_trend.csv"), rd(mini, "stl_trends.csv"), rel_cont)
num_cmp("mini detrended", rd(o, "detrended.csv"), rd(mini, "detrended.csv"), rel_cont)
num_cmp("mini annual maxima", cells_by_year(file.path(o, "annual_maxima.csv.gz")), rd(mini, "annual_maxima.csv"), rel_cont)
g <- rd(o, "gev_params.csv"); w <- rd(mini, "gev_params.csv")
for (p in c("loc", "scale", "shape")) num_cmp(paste("mini GEV", p), g[[p]], w[[p]], rel_cont)
num_cmp("mini frechet", cells_by_year(file.path(o, "frechet.csv.gz")), rd(mini, "frechet.csv"), rel_cont)
num_cmp("mini madogram matrix", rd(o, "madogram_matrix.csv"), rd(mini, "madogram_matrix.csv"), rel15)
s <- rd(o, "edc_summary.csv"); w <- rd(mini, "edc_clustering_params.csv")
num_cmp("mini q30", s$q_threshold[s$order == "r"], w$q30, rel15)
exact_cmp("mini k_edc (R order)", s$k[s$order == "r"], w$k_edc)
exact_cmp("mini k_edc (py order)", s$k[s$order == "py"], w$k_edc)
lab <- rd(o, "labels.csv"); w <- rd(mini, "clusters_edc.csv")
exact_cmp("mini EDC labels (R order)", lab$label_r[order(lab$r_index)], w$cluster)
check("mini EDC labels (py order) same partition", ari(lab$label_py, w$cluster) == 1,
      sprintf("ARI %.17g", ari(lab$label_py, w$cluster)))
# The uint16 condensed file plus the LUT reproduce the madogram exactly.
mm_tok <- tokens(file.path(o, "madogram_matrix.csv"))
gc_ <- rd(o, "grid_cells.csv"); r_of_py <- order(gc_$py_index)
con <- gzfile(file.path(o, "madogram_sum_condensed_py.u16le.bin.gz"), "rb")
S <- readBin(con, "integer", n = 301, size = 2, signed = FALSE, endian = "little"); close(con)
lut <- tokens(file.path(o, "madogram_value_lut.csv"))
lutS <- round(as.numeric(lut[, "S"]))
mpy <- mm_tok[r_of_py, r_of_py]
tok_cmp("mini condensed S + LUT == madogram (py order)", lut[match(S, lutS), "value"], mpy[lower.tri(mpy)])
st <- rd(o, "run_info.csv"); st <- setNames(st$value, st$parameter)
note("mini diagnostics: ties in maxima=", st[["n_cells_with_ties_in_annual_maxima"]],
     ", ties in frechet=", st[["n_cells_with_ties_in_frechet"]],
     ", clipped by R floor=", st[["n_frechet_values_clipped_by_r_floor"]])

# Same chain on the 15-digit CSV, to size the effect of the rounded input.
o2 <- file.path(scratch, "mini_edc_csv")
run("r_edc_fullgrid.R", c(fn, "--monthly-csv", file.path(mini, "monthly_pr.csv"), "--nrow", "5", "--ncol", "5",
                          "--out", o2, "--save-matrix", "--min-k", "2"))
d_tr <- max(abs(as.matrix(rd(o2, "stl_trend.csv")) - as.matrix(rd(mini, "stl_trends.csv"))))
d_fr <- max(abs(cells_by_year(file.path(o2, "frechet.csv.gz")) - as.matrix(rd(mini, "frechet.csv"))))
s2 <- rd(o2, "edc_summary.csv")
note(sprintf("from the 15-digit CSV input instead: trend diff %.2g, frechet diff %.2g, k=%d (info only)",
             d_tr, d_fr, s2$k[s2$order == "r"]))
tok_cmp("mini madogram identical from CSV and exact input",
        tokens(file.path(o2, "madogram_matrix.csv")), mm_tok)

# ── 2. Top-level Saunders/EDC (10x10, n_sim = 10, many tied values) ────────
cat("\n== 2. Saunders EDC (tests/reference_data, %.17e) ==\n")
o <- file.path(scratch, "saunders")
run("r_edc_fullgrid.R", c(fn, "--sim-bin", file.path(ref, "simulation_nonstat_seed42.bin"), "--n-sim", "10",
                          "--nrow", "10", "--ncol", "10", "--out", o, "--save-matrix", "--min-k", "5"))
tok_cmp("saunders madogram matrix (bit-identical)", tokens(file.path(o, "madogram_matrix.csv")),
        tokens(file.path(ref, "madogram_matrix.csv")))
s <- tokens(file.path(o, "edc_summary.csv"), c("q_threshold", "k"))
w <- tokens(file.path(ref, "saunders_clustering_params.csv"))
tok_cmp("saunders q30 (bit-identical)", s[1, "q_threshold"], w[, "q30_threshold"])
tok_cmp("saunders k, exact madogram (R order)", s[1, "k"], w[, "k_saunders"])
# Not a check.  The frozen k comes from R order on the read-back matrix.
# With exact values and the cells permuted, tied merges reorder.
sp <- rd(o, "edc_summary.csv")
note(sprintf("exact madogram, py cell order: k=%d, heights equal to q30: %d (tie sensitivity, info only)",
             sp$k[sp$order == "py"], sp$n_heights_equal_threshold[sp$order == "py"]))
m_exact <- tokens(file.path(o, "hclust_r_merge.csv"))
w_merge <- tokens(file.path(ref, "hclust_saunders_details.csv"))
lab <- rd(o, "labels.csv"); w_lab <- rd(ref, "clusters_saunders.csv")$cluster
note(sprintf("exact madogram vs frozen hclust: %d/99 heights differ (max %.2g), ARI R order %.4f, py order %.4f",
             sum(m_exact[, "height"] != w_merge[, "height"]),
             max(abs(as.numeric(m_exact[, "height"]) - as.numeric(w_merge[, "height"]))),
             ari(lab$label_r[order(lab$r_index)], w_lab), ari(lab$label_py, w_lab)))
# The frozen generator clustered the matrix it read back with read.csv.
# Replaying that read reproduces its hclust input.
o <- file.path(scratch, "saunders_replay")
run("r_edc_fullgrid.R", c(fn, "--madogram-csv", file.path(ref, "madogram_matrix.csv"),
                          "--nrow", "10", "--ncol", "10", "--out", o, "--orders", "r", "--min-k", "5"))
# With an exact parser the replay reads the exact matrix, so its target is
# the exact-input run; with the generator's parser it is the frozen hclust.
o_ex <- file.path(scratch, "saunders")
if (parser_exact) {
  tgt <- "exact-input run"; t_merge <- m_exact; t_order <- tokens(file.path(o_ex, "hclust_r_order.csv"))
  t_lab <- { l <- rd(o_ex, "labels.csv"); l$label_r[order(l$r_index)] }
  t_q <- tokens(file.path(o_ex, "edc_summary.csv"), c("q_threshold", "k"))[1, ]
} else {
  tgt <- "frozen"; t_merge <- w_merge; t_order <- tokens(file.path(ref, "hclust_saunders_order.csv"))
  t_lab <- w_lab; t_q <- c(q_threshold = unname(w[1, "q30_threshold"]), k = unname(w[1, "k_saunders"]))
}
m <- tokens(file.path(o, "hclust_r_merge.csv"))
tok_cmp(paste("saunders replay hclust heights == ", tgt), m[, "height"], t_merge[, "height"])
num_cmp(paste("saunders replay hclust heights (1e-12) ==", tgt), as.numeric(m[, "height"]), as.numeric(t_merge[, "height"]), 1e-12)
tok_cmp(paste("saunders replay hclust merges ==", tgt), m[, c("merge_i", "merge_j")], t_merge[, c("merge_i", "merge_j")])
tok_cmp(paste("saunders replay hclust order ==", tgt), tokens(file.path(o, "hclust_r_order.csv")), t_order)
s <- tokens(file.path(o, "edc_summary.csv"), c("q_threshold", "k"))
tok_cmp(paste("saunders replay q30 ==", tgt), s[1, "q_threshold"], t_q[["q_threshold"]])
tok_cmp(paste("saunders replay k ==", tgt), s[1, "k"], t_q[["k"]])
lab <- rd(o, "labels.csv")
exact_cmp(paste("saunders replay labels ==", tgt), lab$label_r[order(lab$r_index)], t_lab)

# ── 3. Mini CMIP6 LEC on frozen smoothed estimates (15 digits) ─────────────
cat("\n== 3. Mini CMIP6 LEC (cmip6_mini, 15 digits) ==\n")
o <- file.path(scratch, "mini_lec")
run("r_lec_on_estimates.R", c(fn, "--estimates", file.path(mini, "local_estimates_smoothed.csv"),
                              "--out", o, "--quantile", "0.3", "--min-k", "2", "--save-matrix"))
num_cmp("mini ellipse matrix (raw R)", rd(o, "ellipse_raw_matrix.csv"),
        rd(mini, "ellipse_dissimilarity_matrix.csv"), rel15)
s <- rd(o, "lec_summary.csv"); w <- rd(mini, "lec_clustering_params.csv")
num_cmp("mini q30_lec (raw R)", s$q_threshold[s$variant == "raw"], w$q30_lec, rel15)
exact_cmp("mini k_lec (raw R)", s$k[s$variant == "raw"], w$k_lec)
lab <- rd(o, "labels.csv")
exact_cmp("mini LEC labels (raw R)", lab$label_raw, rd(mini, "clusters_lec.csv")$cluster)
note(sprintf("mini corrected variant: q=%.6g k_raw=%d k=%d, empty pairs=%d (info only)",
             s$q_threshold[s$variant == "corrected"], s$k_raw[s$variant == "corrected"],
             s$k[s$variant == "corrected"], s$n_empty_pairs[s$variant == "corrected"]))

# ── 4. Top-level LEC (10x10, %.17e, q10, min k 5) ──────────────────────────
cat("\n== 4. Top-level LEC (tests/reference_data, q10, min k 5) ==\n")
est_bin <- file.path(scratch, "lec10_estimates.f64le.bin")
py <- paste(
  "import csv,struct,sys",
  "r=list(csv.DictReader(open(sys.argv[1])))",
  "v=[float(x[c]) for c in ('a_sm','b_sm','g_sm') for x in r]",
  "open(sys.argv[2],'wb').write(struct.pack('<%dd'%len(v),*v))",
  "print(len(r))", sep = "\n")
n_est <- system2("python3", c("-I", "-c", shQuote(py), file.path(ref, "local_estimates_smoothed.csv"), est_bin), stdout = TRUE)
o <- file.path(scratch, "lec10")
run("r_lec_on_estimates.R", c(fn, "--estimates-bin", est_bin, "--n", n_est, "--out", o,
                              "--quantile", "0.1", "--min-k", "5", "--save-matrix"))
tok_cmp("lec10 ellipse matrix raw (bit-identical)", tokens(file.path(o, "ellipse_raw_matrix.csv")),
        tokens(file.path(ref, "ellipse_dissimilarity_matrix.csv")))
s <- tokens(file.path(o, "lec_summary.csv"), c("q_threshold", "k"))
w <- tokens(file.path(ref, "lec_clustering_params.csv"))
# The frozen q10 was taken on the read-back matrix (generator line 745), so it
# is checked on the replay below.  Here only the size of the difference.
note("lec10 q10 from the exact matrix ", s[1, "q_threshold"], ", frozen ", w[, "q10_threshold"], " (info only)")
tok_cmp("lec10 k_lec, exact matrix", s[1, "k"], w[, "k_lec"])
w_merge <- tokens(file.path(ref, "hclust_lec_details.csv"))
m_exact <- tokens(file.path(o, "hclust_raw_merge.csv"))
note(sprintf("exact matrix vs frozen hclust: %d/99 heights differ as strings (max %.2g), merges differ %d, labels ARI %.4f",
             sum(m_exact[, "height"] != w_merge[, "height"]),
             max(abs(as.numeric(m_exact[, "height"]) - as.numeric(w_merge[, "height"]))),
             sum(m_exact[, c("merge_i", "merge_j")] != w_merge[, c("merge_i", "merge_j")]),
             ari(rd(o, "labels.csv")$label_raw, rd(ref, "clusters_lec.csv")$cluster)))
o <- file.path(scratch, "lec10_replay")
run("r_lec_on_estimates.R", c(fn, "--matrix-csv", file.path(ref, "ellipse_dissimilarity_matrix.csv"),
                              "--out", o, "--quantile", "0.1", "--min-k", "5"))
o_ex <- file.path(scratch, "lec10")
if (parser_exact) {
  tgt <- "exact-input run"; t_merge <- m_exact; t_order <- tokens(file.path(o_ex, "hclust_raw_order.csv"))
  t_lab <- rd(o_ex, "labels.csv")$label_raw
  t_q <- tokens(file.path(o_ex, "lec_summary.csv"), c("q_threshold", "k"))[1, ]
} else {
  tgt <- "frozen"; t_merge <- w_merge; t_order <- tokens(file.path(ref, "hclust_lec_order.csv"))
  t_lab <- rd(ref, "clusters_lec.csv")$cluster; t_q <- c(q_threshold = unname(w[1, "q10_threshold"]), k = unname(w[1, "k_lec"]))
}
m <- tokens(file.path(o, "hclust_given_merge.csv"))
tok_cmp(paste("lec10 replay hclust heights ==", tgt), m[, "height"], t_merge[, "height"])
num_cmp(paste("lec10 replay hclust heights (1e-12) ==", tgt), as.numeric(m[, "height"]), as.numeric(t_merge[, "height"]), 1e-12)
tok_cmp(paste("lec10 replay hclust merges ==", tgt), m[, c("merge_i", "merge_j")], t_merge[, c("merge_i", "merge_j")])
tok_cmp(paste("lec10 replay hclust order ==", tgt), tokens(file.path(o, "hclust_given_order.csv")), t_order)
s <- tokens(file.path(o, "lec_summary.csv"), c("q_threshold", "k"))
tok_cmp(paste("lec10 replay q10 ==", tgt), s[1, "q_threshold"], t_q[["q_threshold"]])
tok_cmp(paste("lec10 replay k_lec ==", tgt), s[1, "k"], t_q[["k"]])
exact_cmp(paste("lec10 replay labels ==", tgt), rd(o, "labels.csv")$label_given, t_lab)

# ── 5. F2 core: R fixed ellipse matrix vs weatherisk, bit for bit ─────────
# Needs the project venv.  Python runs with R hidden (PATH=/usr/bin:/bin).
venv_py <- ".venv/bin/python"
if (file.exists(venv_py)) {
  cat("\n== 5. R fixed calc_distance_ellipses vs weatherisk.clustering ==\n")
  py <- paste(
    "import sys, numpy as np, pandas as pd",
    "sys.path.insert(0, '.')",
    "from weatherisk.clustering import calc_distance_ellipses",
    "e = pd.read_csv(sys.argv[1])[['a_sm','b_sm','g_sm']].to_numpy(np.float64)",
    "e.T.astype('<f8').tofile(sys.argv[2] + '_est.bin')",
    "np.asarray(calc_distance_ellipses(e, res=21), '<f8').tofile(sys.argv[2] + '_py_D.bin')",
    "print(e.shape[0])", sep = "\n")
  for (case in list(c("mini", file.path(mini, "local_estimates_smoothed.csv")),
                    c("lec10", file.path(ref, "local_estimates_smoothed.csv")))) {
    stem <- file.path(scratch, paste0("f2_", case[1]))
    n <- as.integer(system2("env", c("PATH=/usr/bin:/bin", venv_py, "-I", "-c", shQuote(py), case[2], stem),
                            stdout = TRUE))
    o <- paste0(stem, "_r")
    run("r_lec_on_estimates.R", c(fn, "--estimates-bin", paste0(stem, "_est.bin"), "--n", n,
                                  "--out", o, "--variant", "both", "--save-matrix"))
    D_py <- readBin(paste0(stem, "_py_D.bin"), "double", n * n + 1, size = 8, endian = "little")
    tok_cmp(sprintf("%s R fixed ellipse matrix == weatherisk (bit-identical)", case[1]),
            tokens(file.path(o, "ellipse_corrected_matrix.csv")), sprintf("%.17e", D_py))
    s <- rd(o, "lec_summary.csv")
    note(sprintf("%s: empty pairs %d, entries raw vs fixed differing %d", case[1],
                 s$n_empty_pairs[s$variant == "corrected"], s$n_entries_differing_from_raw[s$variant == "corrected"]))
    if (case[1] == "mini")
      exact_cmp("mini raw vs fixed entries differing (44 per tests/test_cmip6_rparity.py)",
                s$n_entries_differing_from_raw[s$variant == "corrected"], 44)
  }
} else note("section 5 skipped: no .venv/bin/python")

res <- do.call(rbind, results)
write.csv(res, file.path(scratch, "validation_results.csv"), row.names = FALSE)
cat(sprintf("\n%d/%d checks passed (%s, %s)\n", sum(res$pass), nrow(res), R.version.string, R.version$arch))
quit(status = if (all(res$pass)) 0 else 1)

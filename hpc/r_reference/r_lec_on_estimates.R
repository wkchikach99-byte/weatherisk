#!/usr/bin/env Rscript
###############################################################################
# R reference for criterion F2: R's LEC clustering stage on given estimates.
#
# Reads smoothed local estimates (a, b, g per cell) and runs, in R:
#   ellipse dissimilarity -> clustering() (hclust average)
#   -> empirical quantile of the upper triangle
#   -> cluster_number_threshold_method -> max(k, min_k) -> cutree.
#
# Two dissimilarity variants:
#   raw        calc_distance_ellipses from r_code/functions.R, unchanged.
#              When both ellipses are empty it runs dist_matrix[j] = 1, a
#              linear index into column 1, so pair [i, j] keeps 0 and cells
#              [1, j] and [j, 1] gain 100 after symmetrisation.
#   corrected  calc_distance_ellipses_fixed from
#              validation_patch_calc_distance_ellipses_fixed.R: the same
#              function with only that line changed to dist_matrix[i,j]=1,
#              which is what weatherisk computes.  It also counts the i<j
#              pairs that reach the empty-ellipse branch (n_empty_pairs).
#              With both variants, lec_summary also reports how many matrix
#              entries and upper-triangle pairs differ between them.
# The two variants are equal when no pair has two empty ellipses.
#
# Row order of the input is kept, so with weatherisk's estimates the cell
# order (and hence hclust tie-breaking) is weatherisk's order.
#
# Arguments:
#   --functions PATH   r_code/functions.R
#   --estimates CSV    columns a_sm,b_sm,g_sm or a,b,g (any other columns are
#                      carried through, e.g. a cell index)
#   --estimates-bin F  the same three columns as little-endian float64,
#                      column-major (all a, then all b, then all g); needs --n.
#                      Use this when the input must be bit-exact: R's decimal
#                      parser is not correctly rounded on aarch64 macOS.
#   --matrix-csv F     a full dissimilarity matrix read with read.csv; skips
#                      the distance step.  Replays the frozen generator, which
#                      read its ellipse matrix back this way.
#   --out DIR
#   --quantile Q       default 0.3
#   --min-k K          default 2 (weatherisk's max(2, k); c4/c5 have no minimum)
#   --res R            ellipse grid resolution, default 21
#   --variant V        raw, corrected or both (default both)
#   --save-matrix      also write the full matrices as CSV (small n only)
#   --patch PATH       the validation patch (default: next to this script)
#
# Outputs in --out (floats as %.17e):
#   lec_summary.csv              per variant q, k_raw, k, empty pairs, timing
#   labels.csv                   row, carried columns, label_<variant>
#   hclust_<variant>_merge.csv   merge_i, merge_j, height
#   hclust_<variant>_order.csv
#   ellipse_<variant>_matrix.csv only with --save-matrix
###############################################################################

t_start <- proc.time()[["elapsed"]]
parse_args <- function(argv) {
  out <- list(); i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    key <- sub("^--", "", a)
    if (grepl("=", key, fixed = TRUE)) {
      out[[sub("=.*$", "", key)]] <- sub("^[^=]*=", "", key); i <- i + 1
    } else if (i < length(argv) && !startsWith(argv[i + 1], "--")) {
      out[[key]] <- argv[i + 1]; i <- i + 2
    } else {
      out[[key]] <- "TRUE"; i <- i + 1
    }
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(name, default = NULL) if (!is.null(opt[[name]])) opt[[name]] else default

functions_path <- get_opt("functions", "r_code/functions.R")
est_path       <- get_opt("estimates")
out_dir        <- get_opt("out")
quantile_level <- as.numeric(get_opt("quantile", "0.3"))
min_k          <- as.integer(get_opt("min-k", "2"))
res            <- as.integer(get_opt("res", "21"))
variant_opt    <- get_opt("variant", "both")
save_matrix    <- identical(get_opt("save-matrix", "FALSE"), "TRUE")
bin_path <- get_opt("estimates-bin"); mat_path <- get_opt("matrix-csv")
if (sum(!is.null(est_path), !is.null(bin_path), !is.null(mat_path)) != 1 || is.null(out_dir))
  stop("give --out and exactly one of --estimates, --estimates-bin, --matrix-csv")
variants <- if (variant_opt == "both") c("raw", "corrected") else variant_opt
if (!all(variants %in% c("raw", "corrected"))) stop("--variant must be raw, corrected or both")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fmt <- function(x) sprintf("%.17e", x)
logmsg <- function(...) {
  cat(sprintf("[%7.1fs] ", proc.time()[["elapsed"]] - t_start), ..., "\n", sep = ""); flush.console()
}

source(functions_path)

# The fixed copy lives in a separate, clearly named validation patch.
patch_path <- get_opt("patch", file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "validation_patch_calc_distance_ellipses_fixed.R"))
source(patch_path)

if (!is.null(est_path)) {
  df <- read.csv(est_path)
  pick <- function(cands) { hit <- intersect(cands, names(df)); if (length(hit) == 0) stop("missing column among ", paste(cands, collapse = "/")); hit[1] }
  cols <- c(pick(c("a_sm", "a")), pick(c("b_sm", "b")), pick(c("g_sm", "g")))
  est <- as.matrix(df[, cols])
  carried <- df[, setdiff(names(df), cols), drop = FALSE]
} else if (!is.null(bin_path)) {
  n_bin <- as.integer(get_opt("n"))
  con <- file(bin_path, "rb")
  x <- readBin(con, "double", n = 3 * n_bin + 1, size = 8, endian = "little"); close(con)
  if (length(x) != 3 * n_bin) stop(bin_path, " holds ", length(x), " doubles, expected ", 3 * n_bin)
  est <- matrix(x, n_bin, 3); cols <- c("a", "b", "g")
  carried <- data.frame(row.names = seq_len(n_bin))
} else {
  D_given <- unname(as.matrix(read.csv(mat_path)))
  est <- matrix(0, nrow(D_given), 3); cols <- "matrix"
  carried <- data.frame(row.names = seq_len(nrow(D_given)))
  variants <- "given"
}
storage.mode(est) <- "double"
if (any(!is.finite(est))) stop("non-finite estimates")
n <- nrow(est)
logmsg("input: ", n, " cells (", paste(cols, collapse = ","), ")")

labels <- data.frame(row = seq_len(n), carried)
rows <- list()
for (v in variants) {
  t0 <- proc.time()[["elapsed"]]
  D <- switch(v, raw = calc_distance_ellipses(est, res = res),
              corrected = calc_distance_ellipses_fixed(est, res = res),
              given = D_given)
  t_dist <- proc.time()[["elapsed"]] - t0
  if (v == "raw" && "corrected" %in% variants) D_raw_keep <- D
  if (v == "corrected" && exists("D_raw_keep")) {
    dd <- D - D_raw_keep
    n_entries_changed <- sum(dd != 0)
    n_upper_pairs_changed <- sum(dd[upper.tri(dd)] != 0)
    rm(D_raw_keep, dd); invisible(gc())
  }
  logmsg(v, ": distance done")
  if (save_matrix) write.csv(matrix(fmt(D), n), file.path(out_dir, paste0("ellipse_", v, "_matrix.csv")),
                             row.names = FALSE, quote = FALSE)
  q_thr <- unname(quantile(D[upper.tri(D)], quantile_level))
  hc <- clustering(D)
  rm(D); invisible(gc())
  k_raw <- cluster_number_threshold_method(hc, q_thr)
  k <- max(k_raw, min_k)
  labels[[paste0("label_", v)]] <- as.integer(cutree(hc, k))
  write.csv(data.frame(merge_i = hc$merge[, 1], merge_j = hc$merge[, 2], height = fmt(hc$height)),
            file.path(out_dir, paste0("hclust_", v, "_merge.csv")), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(order = hc$order), file.path(out_dir, paste0("hclust_", v, "_order.csv")),
            row.names = FALSE, quote = FALSE)
  rows[[v]] <- data.frame(variant = v, n_cells = n, quantile_level = quantile_level,
                          q_threshold = fmt(q_thr), k_raw = k_raw, k = k, min_k = min_k, res = res,
                          n_empty_pairs = if (v == "corrected") validation_n_empty_pairs else NA,
                          n_entries_differing_from_raw = if (v == "corrected" && exists("n_entries_changed")) n_entries_changed else NA,
                          n_upper_pairs_differing_from_raw = if (v == "corrected" && exists("n_upper_pairs_changed")) n_upper_pairs_changed else NA,
                          n_heights_equal_threshold = sum(hc$height == q_thr),
                          seconds_distance = round(t_dist, 1),
                          seconds_total = round(proc.time()[["elapsed"]] - t0, 1))
  logmsg(v, ": q=", fmt(q_thr), " k_raw=", k_raw, " k=", k)
  rm(hc); invisible(gc())
}
summ <- do.call(rbind, rows)
write.csv(summ, file.path(out_dir, "lec_summary.csv"), row.names = FALSE, quote = FALSE)
write.csv(labels, file.path(out_dir, "labels.csv"), row.names = FALSE, quote = FALSE)
print(summ)
logmsg("done (", R.version.string, ")")

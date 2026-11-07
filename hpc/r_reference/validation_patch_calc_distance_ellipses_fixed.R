###############################################################################
# VALIDATION PATCH.  Not part of Thorsten Dickhaus's code and not used by
# weatherisk.  r_code/ is never edited.
#
# Defines calc_distance_ellipses_fixed(), a copy of calc_distance_ellipses
# from r_code/functions.R (line 488) with one line changed.  In the branch
# where both ellipses are empty, line 515 reads
#     dist_matrix[j]=1
# which is a linear index, so it sets element [j, 1] and leaves pair [i, j]
# at 0.  The patch writes
#     dist_matrix[i,j]=1
# which is what the Python implementation computes.  The patched branch also
# counts how many i < j pairs reach it, in validation_n_empty_pairs.
#
# The copy is taken from the sourced function itself, so every other line is
# r_code/functions.R verbatim.  The substitution is checked to hit exactly
# one line, and stops otherwise.
#
# Use: source("r_code/functions.R"); source(this file)
###############################################################################

if (!exists("calc_distance_ellipses", mode = "function"))
  stop("source r_code/functions.R before the validation patch")

validation_n_empty_pairs <- NA_real_

local({
  src <- deparse(calc_distance_ellipses, control = c("keepInteger", "keepNA"))
  orig <- "dist_matrix[j] = 1"
  hit <- which(trimws(src) == orig)
  if (length(hit) != 1)
    stop("validation patch: expected one line '", orig, "' in calc_distance_ellipses, found ", length(hit))
  src[hit] <- sub(orig, paste("{ dist_matrix[i, j] = 1;",
                              "validation_n_empty_pairs <<- validation_n_empty_pairs + 1 }"),
                  src[hit], fixed = TRUE)
  # Reset the counter on entry.
  body_start <- which(trimws(src) == "{")[1]
  src <- append(src, "validation_n_empty_pairs <<- 0", after = body_start)
  f <- eval(parse(text = src))
  environment(f) <- globalenv()
  assign("calc_distance_ellipses_fixed", f, envir = globalenv())
  assign("validation_patch_diff",
         c(paste("-", orig), paste("+", trimws(src[hit + 1]))), envir = globalenv())
})

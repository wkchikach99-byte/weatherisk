/// Jaccard-like ellipse-overlap dissimilarity matrix.
///
/// Port of `clustering.calc_distance_ellipses` from Python.  The Rust
/// version streams the pairwise comparison without materialising the
/// full 3-D boolean tensor in memory (the main memory bottleneck in
/// the Python version).
use rayon::prelude::*;

/// Compute the Jaccard-like ellipse overlap dissimilarity matrix.
///
/// Returns the condensed upper-triangle vector (length n*(n-1)/2),
/// scaled 0–100, in the same order as `scipy.spatial.distance.squareform`.
///
/// Mirrors `clustering.calc_distance_ellipses` exactly.
pub fn calc_distance_ellipses_condensed(
    estimates: &[(f64, f64, f64)], // (a, b, g) per grid point
    res: usize,
) -> Vec<f64> {
    let n = estimates.len();

    // ── Build evaluation points (half-circle for symmetry) ──
    // Mirrors the Python grid construction exactly:
    //   xs = np.repeat(np.linspace(-1, 1, res), res)
    //   ys = np.tile(np.linspace(-1, 1, res), res)
    //   mask = ((xs**2 + ys**2) <= res**2) & ((ys > 0) | (xs > 0))
    let linspace: Vec<f64> = (0..res)
        .map(|i| -1.0 + 2.0 * (i as f64) / ((res - 1) as f64))
        .collect();

    let mut xs_all = Vec::with_capacity(res * res);
    let mut ys_all = Vec::with_capacity(res * res);
    for &xv in &linspace {
        for &yv in &linspace {
            xs_all.push(xv);
            ys_all.push(yv);
        }
    }

    let res_sq = (res * res) as f64;
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for k in 0..xs_all.len() {
        let xv = xs_all[k];
        let yv = ys_all[k];
        if (xv * xv + yv * yv) <= res_sq && (yv > 0.0 || xv > 0.0) {
            xs.push(xv);
            ys.push(yv);
        }
    }
    let n_pts = xs.len();

    // ── Pre-compute sqrt(Q_i) for each estimate ──
    let mut sq: Vec<Vec<f64>> = Vec::with_capacity(n);
    for est in estimates.iter() {
        let (a_i, b_i, g_i) = *est;
        let sg = g_i.sin();
        let cg = g_i.cos();
        let ap = a_i + b_i;

        if a_i == 0.0 && ap == 0.0 {
            sq.push(vec![f64::INFINITY; n_pts]);
            continue;
        }

        let mut row = Vec::with_capacity(n_pts);
        for p in 0..n_pts {
            let xv = xs[p];
            let yv = ys[p];
            let qf = xv * xv * (sg * sg / (a_i * a_i) + cg * cg / (ap * ap))
                + 2.0 * xv * yv * sg * cg * (-1.0 / (a_i * a_i) + 1.0 / (ap * ap))
                + yv * yv * (cg * cg / (a_i * a_i) + sg * sg / (ap * ap));
            row.push(qf.max(0.0).sqrt());
        }
        sq.push(row);
    }

    // Semi-major for each point
    let ab: Vec<f64> = estimates.iter().map(|(a, b, _)| a + b).collect();

    // ── Pairwise dissimilarity (condensed upper triangle) ──
    // Streams pair-by-pair: no 3-D boolean arrays needed.
    let n_pairs = n * (n - 1) / 2;
    let mut condensed = vec![0.0f64; n_pairs];

    // Use Rayon for parallel pair computation
    condensed
        .par_iter_mut()
        .enumerate()
        .for_each(|(flat_idx, out)| {
            // Recover (i, j) from condensed index
            // Row i, col j (j > i), same order as scipy squareform
            let (i, j) = condensed_to_pair(n, flat_idx);

            let mx = ab[i].max(ab[j]).max(1e-300);
            let thr = 1.0 / mx;

            let sq_i = &sq[i];
            let sq_j = &sq[j];

            let mut inter: u64 = 0;
            let mut union: u64 = 0;
            for p in 0..n_pts {
                let in_i = sq_i[p] < thr;
                let in_j = sq_j[p] < thr;
                if in_i && in_j {
                    inter += 1;
                    union += 1;
                } else if in_i || in_j {
                    union += 1;
                }
            }

            let inter_f = inter as f64 + 0.5;
            let union_f = union as f64 + 0.5;

            *out = if union_f == 0.5 {
                100.0
            } else {
                100.0 * (1.0 - inter_f / union_f)
            };
        });

    condensed
}

/// Convert condensed index to (i, j) pair.
///
/// Same ordering as `scipy.spatial.distance.squareform`:
/// index k maps to row i, column j where j > i.
#[inline]
fn condensed_to_pair(n: usize, k: usize) -> (usize, usize) {
    let nf = n as f64;
    let kf = k as f64;
    // i = n - 2 - floor(sqrt(-8k + 4n(n-1) - 7) / 2 - 0.5)
    let i = (nf - 2.0 - ((-8.0 * kf + 4.0 * nf * (nf - 1.0) - 7.0).sqrt() / 2.0 - 0.5).floor())
        as usize;
    // j = k + i + 1 + (n-i)(n-i-1)/2 - n(n-1)/2
    // Rearranged to avoid usize underflow: add before subtracting.
    let ni = n - i;
    let j = k + i + 1 + ni * (ni - 1) / 2 - n * (n - 1) / 2;
    (i, j)
}

/// Full square dissimilarity matrix (for testing / compatibility).
pub fn calc_distance_ellipses_full(
    estimates: &[(f64, f64, f64)],
    res: usize,
) -> Vec<f64> {
    let n = estimates.len();
    let condensed = calc_distance_ellipses_condensed(estimates, res);
    let mut full = vec![0.0f64; n * n];

    let mut k = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            full[i * n + j] = condensed[k];
            full[j * n + i] = condensed[k];
            k += 1;
        }
    }
    full
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_condensed_to_pair_n4() {
        // n=4: pairs are (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
        assert_eq!(condensed_to_pair(4, 0), (0, 1));
        assert_eq!(condensed_to_pair(4, 1), (0, 2));
        assert_eq!(condensed_to_pair(4, 2), (0, 3));
        assert_eq!(condensed_to_pair(4, 3), (1, 2));
        assert_eq!(condensed_to_pair(4, 4), (1, 3));
        assert_eq!(condensed_to_pair(4, 5), (2, 3));
    }

    #[test]
    fn test_identical_ellipses_zero_distance() {
        let est = vec![(1.0, 0.5, 0.0); 3];
        let full = calc_distance_ellipses_full(&est, 21);
        // Diagonal should be 0
        assert_eq!(full[0], 0.0);
        assert_eq!(full[4], 0.0);
        assert_eq!(full[8], 0.0);
        // Off-diagonal should be 0 (identical ellipses)
        assert!(full[1].abs() < 1e-10, "identical ellipses should have 0 dissimilarity");
    }

    #[test]
    fn test_symmetry() {
        let est = vec![(0.5, 0.3, 0.1), (1.0, 0.5, -0.2), (0.8, 0.2, 0.5)];
        let full = calc_distance_ellipses_full(&est, 21);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (full[i * 3 + j] - full[j * 3 + i]).abs() < 1e-12,
                    "matrix should be symmetric"
                );
            }
        }
    }
}

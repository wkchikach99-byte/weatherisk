/// Grid/index helper functions matching the R indexing convention.

/// R helper: dist_x(x1, x2) = x1 - x2
#[inline]
pub fn dist_x(x1: f64, x2: f64) -> f64 {
    x1 - x2
}

/// R helper: dist_y(y1, y2) = y1 - y2
#[inline]
pub fn dist_y(y1: f64, y2: f64) -> f64 {
    y1 - y2
}

/// Column-major linear index from (row, col), 0-based.
#[inline]
pub fn grid_number(i: usize, j: usize, nrow: usize, ncol: usize) -> Option<usize> {
    if i < nrow && j < ncol {
        Some(j * nrow + i)
    } else {
        None
    }
}

/// Inverse of `grid_number`, 0-based.
#[inline]
pub fn number_grid(n: usize, nrow: usize, ncol: usize) -> Option<(usize, usize)> {
    let n_grid = nrow.checked_mul(ncol)?;
    if n < n_grid {
        Some((n % nrow, n / nrow))
    } else {
        None
    }
}

/// Index of the nearest regular-grid point to (x, y), column-major.
pub fn koord_num(x: f64, y: f64, x_ax: &[f64], y_ax: &[f64]) -> Option<usize> {
    if x_ax.is_empty() || y_ax.is_empty() {
        return None;
    }

    let mut best_col = 0usize;
    let mut best_dx = f64::INFINITY;
    for (j, &xv) in x_ax.iter().enumerate() {
        let dx = (xv - x).abs();
        if dx < best_dx {
            best_dx = dx;
            best_col = j;
        }
    }

    let mut best_row = 0usize;
    let mut best_dy = f64::INFINITY;
    for (i, &yv) in y_ax.iter().enumerate() {
        let dy = (yv - y).abs();
        if dy < best_dy {
            best_dy = dy;
            best_row = i;
        }
    }

    grid_number(best_row, best_col, y_ax.len(), x_ax.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_number_roundtrip() {
        let nrow = 5;
        let ncol = 5;
        for n in 0..(nrow * ncol) {
            let (i, j) = number_grid(n, nrow, ncol).unwrap();
            assert_eq!(grid_number(i, j, nrow, ncol), Some(n));
        }
    }

    #[test]
    fn test_grid_number_column_major() {
        assert_eq!(grid_number(0, 0, 5, 5), Some(0));
        assert_eq!(grid_number(4, 0, 5, 5), Some(4));
        assert_eq!(grid_number(0, 4, 5, 5), Some(20));
        assert_eq!(grid_number(4, 4, 5, 5), Some(24));
    }

    #[test]
    fn test_koord_num_regular_grid() {
        let x_ax = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let y_ax = [2.0, 1.0, 0.0, -1.0, -2.0];
        assert_eq!(koord_num(0.0, 0.0, &x_ax, &y_ax), Some(12));
    }
}
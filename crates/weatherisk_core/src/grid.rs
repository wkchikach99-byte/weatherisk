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

    fn parse_csv_fixture(csv_text: &str) -> Vec<Vec<String>> {
        csv_text
            .lines()
            .skip(1)
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split(',')
                    .map(|field| field.trim().trim_matches('"').to_string())
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn load_axis_fixture(csv_text: &str, column_name: &str) -> Vec<f64> {
        let rows = parse_csv_fixture(csv_text);
        let header = csv_text
            .lines()
            .next()
            .expect("axis fixture must have a header row");
        let columns: Vec<&str> = header
            .split(',')
            .map(|field| field.trim().trim_matches('"'))
            .collect();
        let index = columns
            .iter()
            .position(|&name| name == column_name)
            .expect("axis fixture column missing");
        rows.into_iter()
            .map(|row| row[index].parse::<f64>().expect("invalid axis value"))
            .collect()
    }

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

    #[test]
    fn test_dist_helpers_match_r_fixtures() {
        let fixture = include_str!("../../../tests/reference_data/dist_helper_test_cases.csv");

        for row in parse_csv_fixture(fixture) {
            let x1 = row[0].parse::<f64>().expect("invalid x1");
            let x2 = row[1].parse::<f64>().expect("invalid x2");
            let y1 = row[2].parse::<f64>().expect("invalid y1");
            let y2 = row[3].parse::<f64>().expect("invalid y2");
            let expected_dx = row[4].parse::<f64>().expect("invalid dist_x");
            let expected_dy = row[5].parse::<f64>().expect("invalid dist_y");

            assert!((dist_x(x1, x2) - expected_dx).abs() <= 1e-14);
            assert!((dist_y(y1, y2) - expected_dy).abs() <= 1e-14);
        }
    }

    #[test]
    fn test_grid_number_match_r_fixtures() {
        let fixture = include_str!("../../../tests/reference_data/grid_number_test_cases.csv");

        for row in parse_csv_fixture(fixture) {
            let i = row[0].parse::<usize>().expect("invalid i") - 1;
            let j = row[1].parse::<usize>().expect("invalid j") - 1;
            let expected_n = row[2].parse::<usize>().expect("invalid grid_num") - 1;

            assert_eq!(grid_number(i, j, 10, 10), Some(expected_n));
        }
    }

    #[test]
    fn test_number_grid_match_r_fixtures() {
        let fixture = include_str!("../../../tests/reference_data/number_grid_test_cases.csv");

        for row in parse_csv_fixture(fixture) {
            let n = row[0].parse::<usize>().expect("invalid n") - 1;
            let expected_i = row[1].parse::<usize>().expect("invalid i") - 1;
            let expected_j = row[2].parse::<usize>().expect("invalid j") - 1;

            assert_eq!(number_grid(n, 10, 10), Some((expected_i, expected_j)));
        }
    }

    #[test]
    fn test_koord_num_match_r_fixtures() {
        let x_ax = load_axis_fixture(
            include_str!("../../../tests/reference_data/x_axis.csv"),
            "x_ax",
        );
        let y_ax = load_axis_fixture(
            include_str!("../../../tests/reference_data/y_axis.csv"),
            "y_ax",
        );
        let fixture = include_str!("../../../tests/reference_data/koord_num_test_cases.csv");

        for row in parse_csv_fixture(fixture) {
            let x = row[0].parse::<f64>().expect("invalid x");
            let y = row[1].parse::<f64>().expect("invalid y");
            let expected_n = row[2].parse::<usize>().expect("invalid grid_num") - 1;

            assert_eq!(koord_num(x, y, &x_ax, &y_ax), Some(expected_n));
        }
    }
}
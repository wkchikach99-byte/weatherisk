/// weatherisk_core — Rust acceleration for weatherisk.
///
/// Exposes the hot numerical kernels to Python via PyO3/numpy:
///
/// - `neg_log_likelihood_sum`: vectorised objective for L-BFGS-B
/// - `pairwise_density_summand_vec`: element-wise density (for testing)
/// - `cov_fkt_2d`: scalar covariance function (for testing)
/// - `calc_distance_ellipses`: full LEC dissimilarity matrix
/// - `calc_distance_ellipses_condensed`: condensed upper-triangle LEC
/// - `optimize_pairwise_density`: full optimizer loop (global MLE)
/// - `optimize_local_mle`: full optimizer loop (local MLE)
mod density;
mod grid;
mod lec;
mod optimizer;

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

#[pyfunction]
fn dist_x(x1: f64, x2: f64) -> f64 {
    grid::dist_x(x1, x2)
}

#[pyfunction]
fn dist_y(y1: f64, y2: f64) -> f64 {
    grid::dist_y(y1, y2)
}

#[pyfunction]
fn grid_number(i: usize, j: usize, nrow: usize, ncol: usize) -> PyResult<usize> {
    grid::grid_number(i, j, nrow, ncol).ok_or_else(|| {
        PyIndexError::new_err(format!(
            "Index ({i}, {j}) out of bounds for grid {nrow}x{ncol}"
        ))
    })
}

#[pyfunction]
fn number_grid(n: usize, nrow: usize, ncol: usize) -> PyResult<(usize, usize)> {
    grid::number_grid(n, nrow, ncol).ok_or_else(|| {
        PyIndexError::new_err(format!(
            "Index {n} out of bounds for grid size {}",
            nrow * ncol
        ))
    })
}

#[pyfunction]
fn koord_num(
    x: f64,
    y: f64,
    x_ax: PyReadonlyArray1<f64>,
    y_ax: PyReadonlyArray1<f64>,
) -> PyResult<usize> {
    grid::koord_num(x, y, x_ax.as_slice().unwrap(), y_ax.as_slice().unwrap())
        .ok_or_else(|| PyIndexError::new_err("Cannot locate coordinate on empty grid"))
}

/// Scalar covariance function (for parity testing).
#[pyfunction]
fn cov_fkt_2d_scalar(x: f64, y: f64, alpha: f64, a: f64, b: f64, g: f64) -> f64 {
    density::cov_fkt_2d(x, y, alpha, a, b, g)
}

/// Negative log-likelihood sum over paired observation arrays.
///
/// This is the objective function called by the Python L-BFGS-B optimizer.
/// Doing the sum in Rust avoids crossing the Python boundary per element.
#[pyfunction]
fn neg_log_likelihood_sum(
    z1: PyReadonlyArray1<f64>,
    z2: PyReadonlyArray1<f64>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
    df: f64,
    alpha: f64,
    a: f64,
    b: f64,
    g: f64,
) -> f64 {
    density::neg_log_likelihood_sum(
        z1.as_slice().unwrap(),
        z2.as_slice().unwrap(),
        x.as_slice().unwrap(),
        y.as_slice().unwrap(),
        df,
        alpha,
        a,
        b,
        g,
    )
}

/// Element-wise pairwise density summand (for parity testing).
///
/// Returns an array of the same length as the inputs, matching the
/// Python `pairwise_density_summand(z1, z2, x, y, ...)` output.
#[pyfunction]
fn pairwise_density_summand_vec<'py>(
    py: Python<'py>,
    z1: PyReadonlyArray1<f64>,
    z2: PyReadonlyArray1<f64>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
    df: f64,
    alpha: f64,
    a: f64,
    b: f64,
    g: f64,
) -> Bound<'py, PyArray1<f64>> {
    let z1s = z1.as_slice().unwrap();
    let z2s = z2.as_slice().unwrap();
    let xs = x.as_slice().unwrap();
    let ys = y.as_slice().unwrap();
    let n = z1s.len();

    let mut result = Array1::<f64>::zeros(n);
    for i in 0..n {
        result[i] =
            density::pairwise_density_summand_scalar(z1s[i], z2s[i], xs[i], ys[i], df, alpha, a, b, g);
    }
    result.into_pyarray(py)
}

/// Full square LEC dissimilarity matrix.
///
/// Mirrors `clustering.calc_distance_ellipses(estimates, res)`.
/// Returns ndarray shape (n, n), scaled 0–100.
#[pyfunction]
fn calc_distance_ellipses<'py>(
    py: Python<'py>,
    estimates: PyReadonlyArray2<f64>,
    res: usize,
) -> Bound<'py, PyArray2<f64>> {
    let est = estimates.as_array();
    let n = est.nrows();

    let tuples: Vec<(f64, f64, f64)> = (0..n)
        .map(|i| (est[[i, 0]], est[[i, 1]], est[[i, 2]]))
        .collect();

    let flat = lec::calc_distance_ellipses_full(&tuples, res);

    // Reshape flat Vec into 2-D numpy array
    let arr = numpy::ndarray::Array2::from_shape_vec((n, n), flat).unwrap();
    arr.into_pyarray(py)
}

/// Condensed upper-triangle LEC dissimilarity vector.
///
/// Returns 1-D array of length n*(n-1)/2, same order as
/// `scipy.spatial.distance.squareform`.
#[pyfunction]
fn calc_distance_ellipses_condensed<'py>(
    py: Python<'py>,
    estimates: PyReadonlyArray2<f64>,
    res: usize,
) -> Bound<'py, PyArray1<f64>> {
    let est = estimates.as_array();
    let n = est.nrows();

    let tuples: Vec<(f64, f64, f64)> = (0..n)
        .map(|i| (est[[i, 0]], est[[i, 1]], est[[i, 2]]))
        .collect();

    let condensed = lec::calc_distance_ellipses_condensed(&tuples, res);
    Array1::from_vec(condensed).into_pyarray(py)
}

/// The Python module definition.
#[pymodule]
fn weatherisk_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dist_x, m)?)?;
    m.add_function(wrap_pyfunction!(dist_y, m)?)?;
    m.add_function(wrap_pyfunction!(grid_number, m)?)?;
    m.add_function(wrap_pyfunction!(number_grid, m)?)?;
    m.add_function(wrap_pyfunction!(koord_num, m)?)?;
    m.add_function(wrap_pyfunction!(cov_fkt_2d_scalar, m)?)?;
    m.add_function(wrap_pyfunction!(neg_log_likelihood_sum, m)?)?;
    m.add_function(wrap_pyfunction!(nll_with_gradient_py, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_density_summand_vec, m)?)?;
    m.add_function(wrap_pyfunction!(calc_distance_ellipses, m)?)?;
    m.add_function(wrap_pyfunction!(calc_distance_ellipses_condensed, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_pairwise_density_py, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_local_mle_py, m)?)?;
    Ok(())
}

// ── NLL + gradient in one call ───────────────────────────────────────────

/// NLL value and forward-difference gradient in a single FFI crossing.
///
/// Returns (f, numpy array of shape (3,)) so SciPy can use jac=True.
#[pyfunction]
#[pyo3(name = "nll_with_gradient")]
fn nll_with_gradient_py<'py>(
    py: Python<'py>,
    z1: PyReadonlyArray1<f64>,
    z2: PyReadonlyArray1<f64>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
    df: f64,
    alpha: f64,
    a: f64,
    b: f64,
    g: f64,
) -> (f64, Bound<'py, PyArray1<f64>>) {
    let (fval, grad) = density::nll_with_gradient(
        z1.as_slice().unwrap(),
        z2.as_slice().unwrap(),
        x.as_slice().unwrap(),
        y.as_slice().unwrap(),
        df, alpha, a, b, g,
    );
    (fval, Array1::from_vec(grad.to_vec()).into_pyarray(py))
}

// ── Optimizer bindings ───────────────────────────────────────────────────

/// Full multi-start L-BFGS-B for global pairwise density MLE.
///
/// Takes observation matrix z (n_grid × n_sim, row-major), coordinates,
/// and returns optimal [a, b, gamma]. The entire optimizer loop runs in
/// compiled Rust — only one PyO3 boundary crossing.
#[pyfunction]
#[pyo3(name = "optimize_pairwise_density")]
fn optimize_pairwise_density_py<'py>(
    py: Python<'py>,
    z: PyReadonlyArray2<f64>,
    df: f64,
    alpha: f64,
    x_coords: PyReadonlyArray1<f64>,
    y_coords: PyReadonlyArray1<f64>,
    lower_a: f64,
    lower_b: f64,
    upper_a: f64,
    upper_b: f64,
    ensemble: usize,
    max_dist: f64,
    seed: u64,
) -> Bound<'py, PyArray1<f64>> {
    let z_arr = z.as_array();
    let n_grid = z_arr.nrows();
    let n_sim = z_arr.ncols();

    // Flatten to row-major Vec
    let z_flat: Vec<f64> = z_arr.iter().copied().collect();
    let x = x_coords.as_slice().unwrap();
    let y = y_coords.as_slice().unwrap();

    let result = optimizer::optimize_pairwise_density(
        &z_flat, n_grid, n_sim, df, alpha,
        x, y,
        lower_a, lower_b, upper_a, upper_b,
        ensemble, max_dist, seed,
    );

    Array1::from_vec(result.to_vec()).into_pyarray(py)
}

/// Full multi-start L-BFGS-B for local MLE.
///
/// Takes pre-built pair arrays and returns optimal [a, b, gamma].
#[pyfunction]
#[pyo3(name = "optimize_local_mle")]
fn optimize_local_mle_py<'py>(
    py: Python<'py>,
    zi: PyReadonlyArray1<f64>,
    zj: PyReadonlyArray1<f64>,
    xl: PyReadonlyArray1<f64>,
    yl: PyReadonlyArray1<f64>,
    df: f64,
    alpha: f64,
    lower_a: f64,
    lower_b: f64,
    lower_g: f64,
    upper_a: f64,
    upper_b: f64,
    upper_g: f64,
    ensemble: usize,
    seed: u64,
) -> Bound<'py, PyArray1<f64>> {
    let result = optimizer::optimize_local_mle(
        zi.as_slice().unwrap(),
        zj.as_slice().unwrap(),
        xl.as_slice().unwrap(),
        yl.as_slice().unwrap(),
        df, alpha,
        lower_a, lower_b, lower_g,
        upper_a, upper_b, upper_g,
        ensemble, seed,
    );

    Array1::from_vec(result.to_vec()).into_pyarray(py)
}

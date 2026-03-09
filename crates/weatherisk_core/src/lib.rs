/// weatherisk_core — Rust acceleration for weatherisk.
///
/// Exposes the hot numerical kernels to Python via PyO3/numpy:
///
/// - `neg_log_likelihood_sum`: vectorised objective for L-BFGS-B
/// - `pairwise_density_summand_vec`: element-wise density (for testing)
/// - `cov_fkt_2d`: scalar covariance function (for testing)
/// - `calc_distance_ellipses`: full LEC dissimilarity matrix
/// - `calc_distance_ellipses_condensed`: condensed upper-triangle LEC
mod density;
mod lec;

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;

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
    m.add_function(wrap_pyfunction!(cov_fkt_2d_scalar, m)?)?;
    m.add_function(wrap_pyfunction!(neg_log_likelihood_sum, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_density_summand_vec, m)?)?;
    m.add_function(wrap_pyfunction!(calc_distance_ellipses, m)?)?;
    m.add_function(wrap_pyfunction!(calc_distance_ellipses_condensed, m)?)?;
    Ok(())
}

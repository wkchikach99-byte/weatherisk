/// Anisotropic covariance and pairwise composite likelihood density.
///
/// Every expression mirrors `weatherisk/density.py` and
/// `weatherisk/covariance.py` line-for-line to preserve floating-point
/// evaluation order.  Do NOT simplify or factor expressions — the goal
/// is bit-level parity with the Python implementation.
use statrs::function::gamma::ln_gamma;

// ── Covariance ───────────────────────────────────────────────────────────

/// Stationary anisotropic covariance: C(x,y) = exp(-sqrt(Q)^alpha)
///
/// Mirrors `covariance.cov_fkt_2d` exactly.
#[inline(always)]
pub fn cov_fkt_2d(x: f64, y: f64, alpha: f64, a: f64, b: f64, g: f64) -> f64 {
    let sg = g.sin();
    let cg = g.cos();
    let ap = a + b; // semi-major

    let qf = x * x * (sg * sg / (a * a) + cg * cg / (ap * ap))
        + 2.0 * x * y * sg * cg * (-1.0 / (a * a) + 1.0 / (ap * ap))
        + y * y * (cg * cg / (a * a) + sg * sg / (ap * ap));

    (-qf.max(0.0).sqrt().powf(alpha)).exp()
}

// ── Student-t helpers ────────────────────────────────────────────────────

/// Log normalising constant of the Student-t PDF.
///
/// Mirrors `density._t_pdf_log_coeff`.
#[inline(always)]
fn t_pdf_log_coeff(df: f64) -> f64 {
    ln_gamma((df + 1.0) / 2.0) - ln_gamma(df / 2.0) - 0.5 * (df * std::f64::consts::PI).ln()
}

/// Student-t PDF (scalar).
///
/// Mirrors `density._t_pdf`.
#[inline(always)]
fn t_pdf(x: f64, df: f64) -> f64 {
    let coeff = t_pdf_log_coeff(df).exp();
    coeff * (1.0 + (x * x) / df).powf(-(df + 1.0) / 2.0)
}

/// Student-t CDF via the regularised incomplete beta function.
///
/// Mirrors `scipy.special.stdtr(df, x)`.
///
/// Uses the identity:
///   F_t(x; df) = 1 − 0.5 * I_{df/(df+x²)}(df/2, 1/2)   for x > 0
///   F_t(x; df) = 0.5 * I_{df/(df+x²)}(df/2, 1/2)         for x < 0
#[inline(always)]
fn t_cdf(x: f64, df: f64) -> f64 {
    use statrs::function::beta::beta_reg;
    if x == 0.0 {
        return 0.5;
    }
    let t2 = x * x;
    let ib = beta_reg(df / 2.0, 0.5, df / (df + t2));
    if x > 0.0 {
        1.0 - 0.5 * ib
    } else {
        0.5 * ib
    }
}

/// Derivative of the Student-t PDF.
///
/// Mirrors `density._dtdiff`.
#[inline(always)]
fn dtdiff(x: f64, df: f64) -> f64 {
    let pdf = t_pdf(x, df);
    -((df + 1.0) * x / (df + x * x)) * pdf
}

// ── Core density summand ─────────────────────────────────────────────────

/// Log pairwise density contribution for one observation pair.
///
/// Exact port of `density.pairwise_density_summand` from Python.
/// Every sub-expression preserves the same decomposition and evaluation
/// order as the Python/NumPy source to minimise floating-point divergence.
#[inline]
pub fn pairwise_density_summand_scalar(
    z1: f64,
    z2: f64,
    x: f64,
    y: f64,
    df: f64,
    alpha: f64,
    a: f64,
    b: f64,
    g: f64,
) -> f64 {
    let cv = cov_fkt_2d(x, y, alpha, a, b, g);
    // np.sqrt(1 - cv * cv) / np.sqrt(df + 1)
    let c = (1.0 - cv * cv).sqrt() / (df + 1.0).sqrt();

    let inv_df = 1.0 / df;
    let df1 = df + 1.0;

    // m1 = ((z2 / z1) ** (1.0 / df) - cv) / c
    let m1 = ((z2 / z1).powf(inv_df) - cv) / c;
    // m2 = ((z1 / z2) ** (1.0 / df) - cv) / c
    let m2 = ((z1 / z2).powf(inv_df) - cv) / c;

    let dt_m1 = t_pdf(m1, df1);
    let dt_m2 = t_pdf(m2, df1);
    let pt_m1 = t_cdf(m1, df1);
    let pt_m2 = t_cdf(m2, df1);

    // ── First factor (mirrors Python line-for-line) ──
    // term1_a = -pt_m1 / (z1 * z1)
    let term1_a = -pt_m1 / (z1 * z1);
    // term1_b = -dt_m1 * z2 ** (1.0 / df) * z1 ** (-1.0 / df - 2) / c / df
    let term1_b = -dt_m1 * z2.powf(inv_df) * z1.powf(-inv_df - 2.0) / c / df;
    // term1_c = dt_m2 * z1 ** (1.0 / df - 1) * z2 ** (-1.0 / df - 1) / c / df
    let term1_c = dt_m2 * z1.powf(inv_df - 1.0) * z2.powf(-inv_df - 1.0) / c / df;
    let factor1 = term1_a + term1_b + term1_c;

    // ── Second factor ──
    // term2_a = -pt_m2 / (z2 * z2)
    let term2_a = -pt_m2 / (z2 * z2);
    // term2_b = -dt_m2 * z1 ** (1.0 / df) * z2 ** (-1.0 / df - 2) / c / df
    let term2_b = -dt_m2 * z1.powf(inv_df) * z2.powf(-inv_df - 2.0) / c / df;
    // term2_c = dt_m1 * z2 ** (1.0 / df - 1) * z1 ** (-1.0 / df - 1) / c / df
    let term2_c = dt_m1 * z2.powf(inv_df - 1.0) * z1.powf(-inv_df - 1.0) / c / df;
    let factor2 = term2_a + term2_b + term2_c;

    // ── Second-derivative (cross) term ──
    let dtd_m1 = dtdiff(m1, df1);
    let dtd_m2 = dtdiff(m2, df1);

    // Mirrors Python exactly:
    // cross = (
    //     dt_m1 * z1 ** (-1.0 / df - 2) * z2 ** (1.0 / df - 1)
    //     + dt_m2 * z2 ** (-1.0 / df - 2) * z1 ** (1.0 / df - 1)
    //     + dt_m1 * z1 ** (-1.0 / df - 2) * z2 ** (1.0 / df - 1) / df
    //     + dt_m2 * z2 ** (-1.0 / df - 2) * z1 ** (1.0 / df - 1) / df
    //     + dtd_m1 * z1 ** (-2.0 / df - 2) * z2 ** (2.0 / df - 1) / c / df
    //     + dtd_m2 * z2 ** (-2.0 / df - 2) * z1 ** (2.0 / df - 1) / c / df
    // ) / c / df
    let cross = (dt_m1 * z1.powf(-inv_df - 2.0) * z2.powf(inv_df - 1.0)
        + dt_m2 * z2.powf(-inv_df - 2.0) * z1.powf(inv_df - 1.0)
        + dt_m1 * z1.powf(-inv_df - 2.0) * z2.powf(inv_df - 1.0) / df
        + dt_m2 * z2.powf(-inv_df - 2.0) * z1.powf(inv_df - 1.0) / df
        + dtd_m1 * z1.powf(-2.0 * inv_df - 2.0) * z2.powf(2.0 * inv_df - 1.0) / c / df
        + dtd_m2 * z2.powf(-2.0 * inv_df - 2.0) * z1.powf(2.0 * inv_df - 1.0) / c / df)
        / c
        / df;

    // V = pt_m1 / z1 + pt_m2 / z2
    let v = pt_m1 / z1 + pt_m2 / z2;

    // np.log(np.maximum(factor1 * factor2 + cross, 1e-300)) - V
    (factor1 * factor2 + cross).max(1e-300).ln() - v
}

/// Vectorised negative log-likelihood sum over parallel arrays.
///
/// Computes `-sum(pairwise_density_summand(z1[i], z2[i], x[i], y[i], ...))`.
/// This is the function the Python optimizer calls as its objective.
pub fn neg_log_likelihood_sum(
    z1: &[f64],
    z2: &[f64],
    x: &[f64],
    y: &[f64],
    df: f64,
    alpha: f64,
    a: f64,
    b: f64,
    g: f64,
) -> f64 {
    let n = z1.len();
    let mut sum = 0.0;
    for i in 0..n {
        sum += pairwise_density_summand_scalar(z1[i], z2[i], x[i], y[i], df, alpha, a, b, g);
    }
    -sum
}

/// NLL value AND forward-difference gradient in a single call.
///
/// Returns `(f, [df/da, df/db, df/dg])`.
/// This eliminates FFI overhead: SciPy calls this once per iteration
/// instead of 4 separate NLL evaluations (1 for f, 3 for approx_fprime).
pub fn nll_with_gradient(
    z1: &[f64],
    z2: &[f64],
    x: &[f64],
    y: &[f64],
    df: f64,
    alpha: f64,
    a: f64,
    b: f64,
    g: f64,
) -> (f64, [f64; 3]) {
    let f0 = neg_log_likelihood_sum(z1, z2, x, y, df, alpha, a, b, g);

    if !f0.is_finite() {
        return (1e20, [0.0, 0.0, 0.0]);
    }

    let params = [a, b, g];
    let eps = 1e-8;
    let mut grad = [0.0; 3];

    for i in 0..3 {
        let h = eps * params[i].abs().max(1.0);
        let mut p = params;
        p[i] += h;
        let fp = neg_log_likelihood_sum(z1, z2, x, y, df, alpha, p[0], p[1], p[2]);
        grad[i] = if fp.is_finite() { (fp - f0) / h } else { 0.0 };
    }

    (f0, grad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cov_at_origin() {
        let c = cov_fkt_2d(0.0, 0.0, 1.0, 1.0, 1.0, 0.0);
        assert!((c - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_cov_decays() {
        let c1 = cov_fkt_2d(0.5, 0.3, 1.0, 1.0, 0.5, 0.2);
        let c2 = cov_fkt_2d(1.0, 0.6, 1.0, 1.0, 0.5, 0.2);
        assert!(c1 > c2, "covariance should decay with distance");
        assert!(c1 < 1.0);
        assert!(c2 > 0.0);
    }

    #[test]
    fn test_t_cdf_symmetry() {
        let p = t_cdf(0.0, 5.0);
        assert!((p - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_t_cdf_monotone() {
        let p1 = t_cdf(-1.0, 5.0);
        let p2 = t_cdf(0.0, 5.0);
        let p3 = t_cdf(1.0, 5.0);
        assert!(p1 < p2);
        assert!(p2 < p3);
        assert!((p1 + p3 - 1.0).abs() < 1e-12, "symmetry around 0");
    }

    #[test]
    fn test_density_finite() {
        let val = pairwise_density_summand_scalar(1.5, 2.0, 0.5, 0.3, 5.0, 1.0, 0.5, 0.5, 0.1);
        assert!(val.is_finite());
    }

    #[test]
    fn test_neg_llh_sum() {
        let z1 = vec![1.5, 2.0, 1.2];
        let z2 = vec![2.0, 1.5, 1.8];
        let x = vec![0.5, 0.3, 0.7];
        let y = vec![0.3, 0.5, 0.2];
        let nll = neg_log_likelihood_sum(&z1, &z2, &x, &y, 5.0, 1.0, 0.5, 0.5, 0.1);
        assert!(nll.is_finite());
    }
}

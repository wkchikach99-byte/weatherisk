/// Anisotropic covariance and pairwise composite likelihood density.
///
/// Performance-critical inner kernel for the MLE optimizer. Each element
/// of the NLL sum calls into `summand_precomp` which is hot-path optimised:
///
///  - Only 2 `powf` calls per element (base powers z1^(1/df), z2^(1/df));
///    all 12+ other fractional powers are derived via mul/div.
///  - Student-t CDF uses exact closed-form polynomial when df/2 is integer
///    (covers the default df=5), eliminating iterative `beta_reg`.
///  - Student-t PDF coefficient is precomputed once per NLL evaluation.
///  - `dtdiff` is inlined — reuses already-computed t_pdf values.
///  - Covariance is cached per spatial pair (elements from the same pair
///    share identical x,y coordinates in the pipeline data layout).
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

    let d = qf.max(0.0).sqrt();
    // Fast path: alpha=1.0 is the pipeline default — skip powf
    if alpha == 1.0 {
        (-d).exp()
    } else {
        (-d.powf(alpha)).exp()
    }
}

// ── Student-t helpers ────────────────────────────────────────────────────

/// Log normalising constant of the Student-t PDF.
#[inline(always)]
fn t_pdf_log_coeff(df: f64) -> f64 {
    ln_gamma((df + 1.0) / 2.0) - ln_gamma(df / 2.0) - 0.5 * (df * std::f64::consts::PI).ln()
}

/// Student-t PDF (scalar) — reference implementation.
#[allow(dead_code)]
#[inline(always)]
fn t_pdf(x: f64, df: f64) -> f64 {
    let coeff = t_pdf_log_coeff(df).exp();
    coeff * (1.0 + (x * x) / df).powf(-(df + 1.0) / 2.0)
}

/// Student-t PDF with precomputed coefficient.
///
/// Avoids 2× ln_gamma + 1× exp per call. Uses integer-power
/// decomposition when (df+1)/2 is a half-integer (covers all odd df,
/// including the default df=5 → exponent=-3.5).
#[inline(always)]
fn t_pdf_fast(x: f64, df: f64, coeff: f64) -> f64 {
    let u = 1.0 + (x * x) / df;
    let half_exp = (df + 1.0) * 0.5;
    let n = half_exp.floor() as u32;
    let has_half = half_exp - n as f64 > 0.25;

    // u^n by repeated squaring (n is small: 2–5 typically)
    let mut u_n = 1.0_f64;
    let mut base = u;
    let mut e = n;
    while e > 0 {
        if e & 1 == 1 {
            u_n *= base;
        }
        base *= base;
        e >>= 1;
    }

    if has_half {
        coeff / (u_n * u.sqrt())
    } else {
        coeff / u_n
    }
}

/// Student-t CDF — closed-form for integer df/2, fallback otherwise.
///
/// For df where df/2 is a positive integer (e.g. df=6 for the default
/// pipeline parameter df_original=5), the regularised incomplete beta
/// I_z(a, 1/2) has an exact finite-sum form:
///
///   I_z(a, 1/2) = 1 − √(1−z) Σ_{j=0}^{a-1} c_j z^j
///
/// This is 5–8× faster than the iterative `beta_reg` from statrs.
#[inline(always)]
fn t_cdf(x: f64, df: f64) -> f64 {
    if x == 0.0 {
        return 0.5;
    }

    let half_df = df * 0.5;
    let a = half_df.round() as u32;

    if (half_df - a as f64).abs() < 1e-12 && a > 0 && a <= 20 {
        // Closed-form path
        let t2 = x * x;
        let denom = df + t2;
        let z = df / denom;
        let sqrt_omz = x.abs() / denom.sqrt(); // = sqrt(1 − z)

        // Polynomial: Σ c_j z^j  where c_j = (0.5)_j / j!
        let mut poly = 1.0;
        let mut c_j = 1.0_f64;
        let mut z_pow = z;
        for j in 1..a as usize {
            c_j *= (j as f64 - 0.5) / j as f64;
            poly += c_j * z_pow;
            z_pow *= z;
        }

        let ib = 1.0 - sqrt_omz * poly;
        if x > 0.0 {
            1.0 - 0.5 * ib
        } else {
            0.5 * ib
        }
    } else {
        // Fallback: iterative beta_reg
        use statrs::function::beta::beta_reg;
        let t2 = x * x;
        let ib = beta_reg(half_df, 0.5, df / (df + t2));
        if x > 0.0 {
            1.0 - 0.5 * ib
        } else {
            0.5 * ib
        }
    }
}

/// Derivative of the Student-t PDF — reference implementation.
#[allow(dead_code)]
#[inline(always)]
fn dtdiff(x: f64, df: f64) -> f64 {
    let pdf = t_pdf(x, df);
    -((df + 1.0) * x / (df + x * x)) * pdf
}

// ── Core density summand ─────────────────────────────────────────────────

/// Inner summand with precomputed covariance and t-distribution constants.
///
/// This is the hot inner loop.  Only TWO `powf` calls remain (the base
/// powers z1^(1/df) and z2^(1/df)); every other fractional power is
/// derived via multiplication and division of these bases.
#[inline(always)]
fn summand_precomp(
    z1: f64,
    z2: f64,
    cv: f64,
    c: f64,
    inv_df: f64,
    df: f64,
    df1: f64,
    df2: f64,
    t_coeff: f64,
) -> f64 {
    // ── Base powers: only 2 powf calls ──────────────────────────────
    let z1_p = z1.powf(inv_df); // z1^(1/df)
    let z2_p = z2.powf(inv_df); // z2^(1/df)

    let z1_inv = 1.0 / z1;
    let z2_inv = 1.0 / z2;
    let z1_p_inv = 1.0 / z1_p; // z1^(-1/df)
    let z2_p_inv = 1.0 / z2_p; // z2^(-1/df)

    // Derived composite powers (all via mul/div, zero powf)
    let z1_inv2 = z1_inv * z1_inv; // z1^(-2)
    let z2_inv2 = z2_inv * z2_inv; // z2^(-2)

    let z1_pn2 = z1_p_inv * z1_inv2; // z1^(-1/df − 2)
    let z2_pn2 = z2_p_inv * z2_inv2; // z2^(-1/df − 2)
    let z1_pp1 = z1_p * z1_inv; // z1^(1/df − 1)
    let z2_pp1 = z2_p * z2_inv; // z2^(1/df − 1)
    let z1_pn1 = z1_p_inv * z1_inv; // z1^(-1/df − 1)
    let z2_pn1 = z2_p_inv * z2_inv; // z2^(-1/df − 1)

    // Double-exponent composites for cross terms
    let z1_p2n2 = z1_p_inv * z1_p_inv * z1_inv2; // z1^(-2/df − 2)
    let z2_p2n2 = z2_p_inv * z2_p_inv * z2_inv2; // z2^(-2/df − 2)
    let z1_p2_inv = z1_p * z1_p * z1_inv; // z1^(2/df − 1)
    let z2_p2_inv = z2_p * z2_p * z2_inv; // z2^(2/df − 1)

    // ── m values (using base powers, no extra powf) ─────────────────
    let m1 = (z2_p * z1_p_inv - cv) / c;
    let m2 = (z1_p * z2_p_inv - cv) / c;

    // ── t-distribution values (precomputed coefficient) ─────────────
    let dt_m1 = t_pdf_fast(m1, df1, t_coeff);
    let dt_m2 = t_pdf_fast(m2, df1, t_coeff);
    let pt_m1 = t_cdf(m1, df1);
    let pt_m2 = t_cdf(m2, df1);

    // ── Inline dtdiff: reuse dt_m1, dt_m2 (avoids 2 redundant t_pdf) ──
    // dtdiff(m, df1) = -((df1+1) * m / (df1 + m²)) * t_pdf(m, df1)
    let dtd_m1 = -(df2 * m1 / (df1 + m1 * m1)) * dt_m1;
    let dtd_m2 = -(df2 * m2 / (df1 + m2 * m2)) * dt_m2;

    // ── factor1 ─────────────────────────────────────────────────────
    let term1_a = -pt_m1 * z1_inv2;
    let term1_b = -dt_m1 * z2_p * z1_pn2 / c / df;
    let term1_c = dt_m2 * z1_pp1 * z2_pn1 / c / df;
    let factor1 = term1_a + term1_b + term1_c;

    // ── factor2 ─────────────────────────────────────────────────────
    let term2_a = -pt_m2 * z2_inv2;
    let term2_b = -dt_m2 * z1_p * z2_pn2 / c / df;
    let term2_c = dt_m1 * z2_pp1 * z1_pn1 / c / df;
    let factor2 = term2_a + term2_b + term2_c;

    // ── Cross term ──────────────────────────────────────────────────
    let cross = (dt_m1 * z1_pn2 * z2_pp1
        + dt_m2 * z2_pn2 * z1_pp1
        + dt_m1 * z1_pn2 * z2_pp1 / df
        + dt_m2 * z2_pn2 * z1_pp1 / df
        + dtd_m1 * z1_p2n2 * z2_p2_inv / c / df
        + dtd_m2 * z2_p2n2 * z1_p2_inv / c / df)
        / c
        / df;

    // V = pt_m1 / z1 + pt_m2 / z2
    let v = pt_m1 * z1_inv + pt_m2 * z2_inv;

    (factor1 * factor2 + cross).max(1e-300).ln() - v
}

/// Log pairwise density contribution for one observation pair.
///
/// Public API for element-wise evaluation and parity tests.
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
    let df1 = df + 1.0;
    let c = (1.0 - cv * cv).sqrt() / df1.sqrt();
    let inv_df = 1.0 / df;
    let df2 = df + 2.0;
    let t_coeff = t_pdf_log_coeff(df1).exp();
    summand_precomp(z1, z2, cv, c, inv_df, df, df1, df2, t_coeff)
}

/// Vectorised negative log-likelihood sum over parallel arrays.
///
/// Computes `-sum(pairwise_density_summand(z1[i], z2[i], x[i], y[i], ...))`.
///
/// **Caching**: Elements from the same spatial pair share identical x,y
/// coordinates (due to `np.repeat` in the Python pair-building code).
/// The covariance is computed once per unique (x,y) block and reused
/// for all temporal observations in that pair — typically a 48× reduction.
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
    let inv_df = 1.0 / df;
    let df1 = df + 1.0;
    let df2 = df + 2.0;
    let t_coeff = t_pdf_log_coeff(df1).exp();
    let sqrt_df1 = df1.sqrt();

    let mut sum = 0.0;
    let mut prev_x = f64::NAN;
    let mut prev_y = f64::NAN;
    let mut cv = 0.0_f64;
    let mut c = 0.0_f64;

    for i in 0..n {
        // Cache covariance per unique (x,y) pair
        if x[i] != prev_x || y[i] != prev_y {
            prev_x = x[i];
            prev_y = y[i];
            cv = cov_fkt_2d(x[i], y[i], alpha, a, b, g);
            c = (1.0 - cv * cv).sqrt() / sqrt_df1;
        }
        sum += summand_precomp(z1[i], z2[i], cv, c, inv_df, df, df1, df2, t_coeff);
    }
    -sum
}

/// NLL value AND forward-difference gradient in a single call.
///
/// Returns `(f, [df/da, df/db, df/dg])`.
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

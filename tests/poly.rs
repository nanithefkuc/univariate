//! The public polynomial surface against the naive oracles.

use fgf::field::{Elem, Field};
use fgf::kernel::FieldKernels;
use fgf::{FanPaar8, FanPaar16, FanPaar32, FanPaar64, Gf8B, Gf16, Gf32, Gf64};
use univariate::{Polynomial, PolynomialError, karatsuba_multiply};

mod oracles;
use oracles::{naive_evaluate, naive_gcd_ext, naive_multiply, noise, noise_poly};

fn assert_ring_identities<F: FieldKernels>() {
    let p = oracles::noise_poly::<F>(11, 0x5EED_0001);
    let q = oracles::noise_poly::<F>(7, 0x5EED_0002);
    let scale = noise::<F>(1, 0x5EED_0003)[0];

    // Product tiers: schoolbook, Karatsuba, and the naive convolution are
    // byte-identical.
    let schoolbook = p
        .multiply_truncated(&q, p.coefficient_count() + q.coefficient_count())
        .expect("schoolbook product");
    let dispatched = p.multiply(&q).expect("dispatched product");
    assert_eq!(schoolbook, naive_multiply(&p, &q));
    assert_eq!(dispatched, naive_multiply(&p, &q));
    let karatsuba = karatsuba_multiply(&p, &q).expect("karatsuba product");
    assert_eq!(karatsuba, naive_multiply(&p, &q));

    // Additive identities at random points.
    let sum = p.add(&q).expect("sum");
    let axpy = p.add_scaled(scale, &q).expect("axpy");
    for seed in 0..5 {
        let point = noise::<F>(1, 0x7000 + seed)[0];
        assert_eq!(
            sum.evaluate(point),
            naive_evaluate(&p, point).add(naive_evaluate(&q, point))
        );
        assert_eq!(
            axpy.evaluate(point),
            naive_evaluate(&p, point).add(scale.mul(naive_evaluate(&q, point)))
        );
        assert_eq!(
            dispatched.evaluate(point),
            naive_evaluate(&p, point).mul(naive_evaluate(&q, point))
        );
        assert_eq!(
            p.shifted(3).expect("shift").evaluate(point),
            point.pow(3).mul(naive_evaluate(&p, point))
        );
    }

    // Division identity: a = q·b + r with deg r < deg b, against the naive
    // multiply for the check.
    let (quotient, remainder) = p.div_rem(&q).expect("division");
    assert!(
        remainder
            .degree()
            .is_none_or(|degree| degree < q.degree().unwrap_or(0))
            || q.is_zero(),
        "remainder degree must fall below the divisor"
    );
    let reconstructed = naive_multiply(&quotient, &q);
    assert_eq!(reconstructed.add(&remainder).expect("reconstruction"), p);

    // Exact division round trip.
    let product = p.multiply(&q).expect("product");
    assert_eq!(product.exact_divide(&p).expect("exact"), q);
    assert_eq!(product.exact_divide(&q).expect("exact"), p);
    assert_eq!(p.exact_divide(&q), Err(PolynomialError::NonExactDivision));

    // gcd cases: g | a, g | b, and the fixed extremes.
    let common = Polynomial::from_coefficients(&[scale, F::Elem::ONE]).expect("common factor");
    let left = common
        .multiply(&oracles::noise_poly::<F>(5, 0x5EED_0010))
        .expect("left");
    let right = common
        .multiply(&oracles::noise_poly::<F>(3, 0x5EED_0011))
        .expect("right");
    let g = left.gcd(&right).expect("gcd");
    assert!(g.degree().is_some_and(|degree| degree >= 1));
    assert!(left.div_rem(&g).expect("div").1.is_zero());
    assert!(right.div_rem(&g).expect("div").1.is_zero());
    assert_eq!(p.gcd(&Polynomial::zero()).expect("gcd"), p.monic());
    assert_eq!(p.gcd(&p).expect("gcd"), p.monic());

    // Extended gcd: Bézout identity via the naive multiply, cofactor bounds,
    // and agreement with the naive extended Euclid.
    let relation = p.gcd_ext(&q).expect("extended gcd");
    let identity = naive_multiply(&relation.a_cofactor, &p)
        .add(&naive_multiply(&relation.b_cofactor, &q))
        .expect("bézout");
    assert_eq!(identity, relation.gcd);
    assert!(relation.gcd.is_zero() || relation.gcd.monic() == relation.gcd);
    let bound = |coefficient: &Polynomial<F>, limit: Option<usize>| {
        coefficient
            .degree()
            .is_none_or(|degree| limit.is_some_and(|limit| degree <= limit))
    };
    let gcd_degree = relation.gcd.degree().unwrap_or(0);
    assert!(bound(
        &relation.a_cofactor,
        q.degree().and_then(|degree| degree.checked_sub(gcd_degree))
    ));
    assert!(bound(
        &relation.b_cofactor,
        p.degree().and_then(|degree| degree.checked_sub(gcd_degree))
    ));
    let (naive_g, naive_s, naive_t) = naive_gcd_ext(&p, &q);
    assert_eq!(
        (relation.gcd, relation.a_cofactor, relation.b_cofactor),
        (naive_g, naive_s, naive_t)
    );

    // Characteristic-two square: bit-spread equals the naive self-product.
    let squared = p.square().expect("square");
    assert_eq!(squared, naive_multiply(&p, &p));

    // Derivatives: char-2 formal derivative of X^n is 0 for even n, X^(n-1)
    // for odd n.
    let mut powers = vec![F::Elem::ZERO; 8];
    powers[7] = F::Elem::ONE;
    let x7 = Polynomial::<F>::from_coefficients(&powers).expect("X^7");
    let derivative = x7.formal_derivative().expect("derivative");
    assert_eq!(derivative.coefficient(6), F::Elem::ONE);
    assert_eq!(derivative.coefficient_count(), 7);

    // Canonical form: no constructor or mutating op leaves a high zero.
    let mut truncated = product.clone();
    truncated.truncate(3);
    assert!(
        truncated.is_zero()
            || truncated.coefficient(truncated.degree().unwrap()).is_one()
            || !truncated.coefficient(truncated.degree().unwrap()).is_zero()
    );
    let zeroed = oracles::noise_poly::<F>(4, 0x5EED_0020).scaled(F::Elem::ZERO);
    assert!(zeroed.is_zero() && zeroed.coefficient_count() == 0);
}

#[test]
fn ring_identities_hold_across_every_field() {
    assert_ring_identities::<Gf8B>();
    assert_ring_identities::<Gf8B>();
    assert_ring_identities::<Gf16>();
    assert_ring_identities::<Gf32>();
    assert_ring_identities::<Gf64>();
    assert_ring_identities::<FanPaar8>();
    assert_ring_identities::<FanPaar16>();
    assert_ring_identities::<FanPaar32>();
    assert_ring_identities::<FanPaar64>();
}

#[test]
fn zero_and_unit_edges() {
    let zero = Polynomial::<Gf16>::zero();
    let one = Polynomial::<Gf16>::one().expect("one");
    let p = noise_poly::<Gf16>(4, 0x5EED_0030);

    assert_eq!(zero.degree(), None);
    assert_eq!(zero.monic(), zero);
    assert_eq!(zero.multiply(&p).expect("product"), zero);
    assert_eq!(p.multiply(&zero).expect("product"), zero);
    assert_eq!(one.multiply(&p).expect("product"), p);
    assert_eq!(p.div_rem(&zero), Err(PolynomialError::DivisionByZero));
    assert_eq!(p.exact_divide(&zero), Err(PolynomialError::DivisionByZero));
    assert_eq!(p.remainder(&one).expect("remainder"), zero);
    // deg a < deg b: quotient zero, remainder a.
    let small = noise_poly::<Gf16>(2, 0x5EED_0031);
    let (quotient, remainder) = small.div_rem(&p).expect("division");
    assert_eq!(quotient, zero);
    assert_eq!(remainder, small);
    // pow_mod with zero modulus is the explicit zero-divisor error (U6).
    assert_eq!(one.pow_mod(1, &zero), Err(PolynomialError::DivisionByZero));
    // inv(0) == 0 is inherited, never inferred as a root or pivot.
    assert!(<Gf16 as Field>::Elem::ZERO.inv().is_zero());
    let zero_leading = Polynomial::<Gf16>::from_coefficients(&[
        <Gf16 as Field>::Elem::ONE,
        <Gf16 as Field>::Elem::ZERO,
        <Gf16 as Field>::Elem::ZERO,
    ])
    .expect("x");
    assert!(zero_leading.coefficient(2).is_zero());
}

#[test]
fn karatsuba_matches_schoolbook_at_every_size_bucket() {
    for len in [16, 17, 31, 47, 48, 49, 64, 96, 128, 200, 300] {
        let left = noise_poly::<Gf8B>(len, 0x5EED_0040);
        let right = noise_poly::<Gf8B>(len.saturating_sub(1).max(1), 0x5EED_0041);
        let schoolbook = naive_multiply(&left, &right);
        assert_eq!(
            karatsuba_multiply(&left, &right).expect("karatsuba"),
            schoolbook,
            "karatsuba mismatch at operand length {len}"
        );
    }
}

#[test]
fn division_error_paths_preserve_state() {
    let dividend = noise_poly::<Gf16>(6, 0x5EED_0050);
    let _divisor = noise_poly::<Gf16>(2, 0x5EED_0051);
    let mut quotient = Polynomial::<Gf16>::one().expect("one");
    let mut remainder = noise_poly::<Gf16>(3, 0x5EED_0052);

    // Division by zero into reusable outputs leaves them untouched.
    let quotient_before = quotient.clone();
    let remainder_before = remainder.clone();
    assert_eq!(
        dividend.div_rem_into(&Polynomial::zero(), &mut quotient, &mut remainder),
        Err(PolynomialError::DivisionByZero)
    );
    assert_eq!(quotient, quotient_before);
    assert_eq!(remainder, remainder_before);

    // A non-exact X-power division leaves the output untouched.
    let value = Polynomial::<Gf16>::from_coefficients(&[
        <Gf16 as Field>::Elem::ONE,
        <Gf16 as Field>::Elem::ZERO,
        <Gf16 as Field>::Elem::ONE,
    ])
    .expect("1 + X^2");
    let mut out = Polynomial::<Gf16>::one().expect("one");
    let out_before = out.clone();
    assert_eq!(
        value.divide_by_x_power_into(1, &mut out),
        Err(PolynomialError::NonExactDivision)
    );
    assert_eq!(out, out_before);
}

#[test]
fn packed_buffer_round_trips() {
    let p = noise_poly::<Gf8B>(9, 0x5EED_0060);
    let packed = p.as_packed().to_vec();
    let rebuilt = Polynomial::<Gf8B>::from_packed(packed).expect("aligned");
    assert_eq!(rebuilt, p);
    // A partial trailing element is rejected: only possible with multi-byte
    // elements, so use a 2-byte field.
    let wide = noise_poly::<Gf16>(5, 0x5EED_0061);
    let mut ragged = wide.as_packed().to_vec();
    ragged.push(0);
    assert!(Polynomial::<Gf16>::from_packed(ragged).is_none());
    // A trailing zero coefficient never survives construction.
    let mut padded = p.as_packed().to_vec();
    padded.extend_from_slice(&[0; 2]);
    let normalized = Polynomial::<Gf8B>::from_packed(padded).expect("aligned");
    assert_eq!(normalized, p);
}

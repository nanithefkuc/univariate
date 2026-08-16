# Changelog

All notable changes to this project are documented in this file.

## 0.0.0 (2026-08-17)

Initial implementation of the univariate polynomial ring over GF(2^m).

- `Polynomial<F>`: dense monomial-basis coefficients packed in `fgf`'s
  little-endian element representation, canonical form enforced everywhere,
  zero polynomial as the empty buffer.
- Ring operations: add / add-scaled / scale / shift, schoolbook product
  through `fgf`'s packed kernels, characteristic-two `O(deg)` squaring,
  Hasse and formal derivatives, affine composition.
- Product tiers: measured schoolbook↔Karatsuba dispatch and the
  `butterfly-fft`-composed AFFT batched product behind the `fft` feature.
- Division: `div_rem` / `exact_divide` / `remainder` / `monic`, `X^k`
  valuation and exact division, `multiply_mod` / `square_mod` / `pow_mod`,
  all with reusable-output forms.
- `gcd`, `gcd_ext` with Bézout cofactors, and `truncated_eea` — the
  key-equation / Padé primitive equivalent to Berlekamp–Massey.
- Truncated power series: Newton-doubling `inverse_mod_x_power`, series
  division, reversal.
- Root finding: classical Chien search, `gcd(p, X^|F|+X)` equal-degree
  extraction with deterministic trace splitting, a standalone
  linearized/affine solver, and Roth–Ruckenstein / Alekhnovich
  power-series lifting over bivariate Y-rows.
- Evaluation and interpolation: Horner, subproduct-tree multipoint over
  arbitrary points, Newton and Lagrange interpolation, `EvaluationDomain`
  over arbitrary / additive-subspace / affine-coset point sets, and
  `butterfly-fft` transform composition under the `fft` feature.
- Measured backend selectors in `cost`; crossovers recorded in
  `BENCHMARKS.md`.

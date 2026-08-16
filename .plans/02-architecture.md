# Architecture

## Module layout

```
src/
  lib.rs            // declarations, re-exports, lint header, scope boundary. No bodies.
  error.rs          // PolynomialError, DomainError — small enums per domain, value+limit

  poly/
    mod.rs
    dense.rs        // Polynomial<F> (packed Vec<u8>, canonical), constructors, accessors
    ring.rs         // add/scale/shift/multiply(schoolbook)/multiply_truncated/square/derivative
    karatsuba.rs    // middle-tier multiply (GAP)
    divide.rs       // div_rem(_into), exact_divide, remainder, monic, x_valuation, *_mod
    gcd.rs          // gcd, extended gcd w/ cofactors, truncated/partial EEA (GAP)
    series.rs       // truncated power series: inverse mod x^t (Newton), division, reciprocal (GAP)

  eval/
    mod.rs
    horner.rs       // evaluate / evaluate_many (single + small point sets)
    multipoint.rs   // subproduct-tree eval + interpolation over ARBITRARY points (GAP for arbitrary)
    newton.rs       // Newton basis + incremental interpolation (from gs-engine)
    domain.rs       // EvaluationDomain<F> + backend selection
    transform.rs    // #[cfg(feature="fft")] compose butterfly-fft: subspace eval/interp, monomial<->novel

  roots/
    mod.rs
    chien.rs        // classical Chien search (GAP)
    linearized.rs   // affine / 2-linearized root solver (GAP)
    equal_degree.rs // base_field_roots: gcd(p, X^|F|+X)+split (from gs-engine)
    lift.rs         // Roth-Ruckenstein / Alekhnovich power-series lifting (from gs-engine)

  cost.rs           // crossover thresholds + backend selectors (from gs-engine)
```

`lib.rs` and every `mod.rs` hold declarations only. No function bodies, no
`impl` blocks. All public items are re-exported at the crate root.

## The `Polynomial<F>` type

One type, packed and canonical (U3). It is the object the whole crate is about;
every operation returns or consumes it.

```rust
pub struct Polynomial<F: FieldKernels> {
    coeffs: Vec<u8>,          // F::BYTES-packed, little-endian, low degree first
    field: PhantomData<F>,    // no per-instance field state; kernels are static
}
```

Three properties, all type invariants and all enforced at construction and after
every mutating op:

1. **Packed, little-endian, low-degree-first.** `coeffs[i*F::BYTES .. (i+1)*
   F::BYTES]` is the coefficient of `x^i`, encoded exactly as `fgf` encodes an
   element. A row is therefore driven straight through `fgf::ops` with no
   transcoding — this is why the storage matches `fgf`'s buffer contract rather
   than a `Vec<F::Elem>` (lifted from `gs-engine/src/poly/mod.rs:56-243`).
2. **Always normalized.** `normalize()` trims trailing zero coefficients after
   every operation that can create them (add, subtract-equals, division). No
   `Polynomial` ever carries a high-order zero.
3. **Zero is the empty buffer.** `coeffs.is_empty()` is the canonical zero
   polynomial. `degree()` of zero is a documented sentinel (`Option<usize>` /
   `usize::MAX`-style), never `0`, so a degree-0 nonzero constant is
   distinguishable from zero. `inv(0) == 0` is inherited from `fgf`
   (`fgf/src/field/mod.rs:45-124`), so leading-coefficient and root tests call
   `is_zero()` explicitly (U6) rather than reading it off a division result.

The API surface splits across `poly/` by concern: `ring.rs` (add/scale/shift/
schoolbook multiply/`multiply_truncated`/char-2 `square_into`/`formal_derivative`/
`hasse_derivative`), `karatsuba.rs` (the GAP middle tier), `divide.rs`
(`div_rem`/`exact_divide`/`remainder`/`monic`/`*_mod`/`x_valuation`), `gcd.rs`
(plain gcd + the GAP cofactor and truncated forms), and `series.rs` (the GAP
truncated-power-series inverse/division/reciprocal).

**Why packed `Vec<u8>`, not `Vec<F::Elem>`.** The whole reason this ring is fast
is that `fgf::ops` processes a run of coefficients as one packed buffer above a
measured lane-bytes crossover (U1, `gs-engine/src/poly/arithmetic.rs:573`
`use_packed_kernel`). A `Vec<F::Elem>` would force a scalar loop or a copy into a
byte buffer before every kernel call. Storing the byte buffer directly makes the
packed kernel the default path and the scalar loop the small-degree exception.

**Scratch, not allocation, on the hot path (U5).** Every op that can allocate has
a `*_into` form writing into a caller-owned `Polynomial` or a per-op `*Scratch`
(`PolynomialProductScratch`, `AlekhnovichScratch`, and the like, lifted from
`gs-engine`). `tests/zero_alloc.rs` proves steady-state zero allocation under a
counting global allocator.

## `EvaluationDomain<F>` — the point-set abstraction

Evaluation and interpolation are parameterized by *where* the points live,
because the point set decides the backend. `EvaluationDomain<F>` is the three
cases the stack actually uses, lifted from `gs-engine/src/domain.rs`:

```rust
pub enum EvaluationDomain<F: FieldKernels> {
    /// Any explicit, ordered point set. Horner or subproduct tree.
    Arbitrary { points: Vec<u8> },              // F::BYTES-packed
    /// A size-2^k additive subspace. Composes butterfly-fft when `fft` is on.
    Subspace { log_size: u32 },
    /// An affine coset alpha + V of an additive subspace.
    AffineCoset { shift: [u8; MAX_ELEM_BYTES], log_size: u32 },
}
```

The domain plus a measured cost model (`cost.rs`) selects the evaluation and
interpolation backend:

- **`Arbitrary`** → Horner for one or a handful of points (`eval/horner.rs`); the
  subproduct-tree multipoint path (`eval/multipoint.rs`) above the measured
  crossover. Interpolation is Newton (`eval/newton.rs`) for small `n`
  (`MODULE_INTERPOLATION_CROSSOVER=8`) and Lagrange-via-subproduct-tree above it.
- **`Subspace` / `AffineCoset`** → `#[cfg(feature="fft")]` composes
  `butterfly-fft`: `monomial_to_novel` then `TransformPlan::forward` for
  evaluation, `TransformPlan::inverse` then `novel_to_monomial` for
  interpolation, `ShiftedPlan` for the coset (`eval/transform.rs`). Without the
  `fft` feature these variants fall back to the arbitrary path over the
  enumerated subspace, so the domain type is total in both build configurations.

Ordering inside a domain is a wire property (U8): the enumeration order of an
`Arbitrary` point set and the index order the transform produces are documented
and frozen, because consumers map roots and interpolated values back to
positions.

## Error model

Hand-rolled, no `thiserror`. One `src/error.rs`, small enums per failure domain,
never a god-enum.

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum PolynomialError {
    /// A construction or operand shape is invalid (e.g. ragged coefficient buffer).
    Config { len: usize, element_bytes: usize },
    /// Division or modular reduction by the zero polynomial.
    DivisionByZero,
    /// exact_divide was asked for a quotient that leaves a nonzero remainder.
    NonExactDivision { remainder_degree: usize },
    /// An operation exceeded a degree/length limit it validates at the boundary.
    Degree { value: usize, limit: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DomainError {
    /// A point set / value vector length does not match the domain size.
    LengthMismatch { expected: usize, found: usize },
    /// A claimed additive subspace / coset is not one.
    NotSubspace { log_size: u32 },
}
```

Every struct variant carries **both the offending value and the limit**. Manual
`impl Display` with inline-captured args; `impl std::error::Error` under `std`.
Every fallible public constructor and operation has a `/// # Errors` section.

Note what is *not* an error: a zero result from `inv` (U6, `inv(0) == 0` is
inherited from `fgf` and is total), a polynomial with no roots, an empty gcd, or
a degree-0 gcd. Only genuine geometry violations (ragged buffers, mismatched
domain lengths), division by the zero polynomial, and a non-exact `exact_divide`
are errors. `PolynomialError::Degree` and `DomainError` are checked at the public
boundary (U9); internal kernels reach the packed `fgf::ops` calls only with
valid geometry and use `debug_assert!`.

Every error path validates before mutation. Tests compare a state snapshot before
and after malformed input, so "returned an error after a partial in-place
division" cannot pass.

## Trait seams

Two bounds, and no new trait of our own.

```rust
// The field bound every polynomial is generic over.
F: fgf::kernel::FieldKernels

// The stronger bound the transform paths require, only where `fft` is on.
#[cfg(feature = "fft")]
F: butterfly_fft::ButterflyKernels    // : FieldKernels + Sealed
```

`FieldKernels` is sealed upstream (`fgf/src/kernel/mod.rs`); `univariate` cannot
add a field, and does not try — new element shapes are composed from the public
`fgf::ops` functions or upstreamed. The transform-carrying functions in
`eval/transform.rs` take the tighter `ButterflyKernels` bound, which
`butterfly-fft` seals to the vectorizable fields (GF(2^8)/GF(2^16); wider fields
fall to scalar upstream, `butterfly-fft` kernel/mod.rs:153). Because that bound
lives only behind `#[cfg(feature="fft")]`, the `--no-default-features` API never
mentions it.

There is **no trait for "a polynomial"** and **no trait for "a field element".**
`Polynomial<F>` is a concrete type; a caller that wants to abstract over rings
writes its own trait, and so far no caller does. The dispatch enums
(`ProductStrategy{Auto, Schoolbook, Afft}` plus the Karatsuba tier,
`RootBackend`, and the interpolation selector) are concrete data read by
`cost.rs`, not trait objects.

## Backend handling

`univariate` owns no kernels and mints no backend policy. `cost.rs` is a
pass-through: where a path needs to know the active backend to pick a blocking
factor or a crossover, it reads `fgf::backend_for::<F>()` and does not re-probe
target features.

Critically, `univariate` introduces **no environment override.** The one
stack-wide, downgrade-only override is `SIMD_BACKEND`, resolved inside
`simdispatch` and re-exported as `fgf::Backend`
(`gfm/.plans/09-fgf-rename.md`). A per-crate `UNIVARIATE_BACKEND` would be a
fourth copy of downgrade-only resolution logic and is exactly the failure `fgf`'s
own planning warns against. There is no such variable.

What `cost.rs` actually decides is the crossover between tiers, and those
thresholds are reused from `gs-engine` rather than re-derived (U7): the AFFT
product crossovers (`AFFT_*_CROSSOVER`), `MODULE_INTERPOLATION_CROSSOVER=8`,
`DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER=20000`, and the product/root/interpolation
selectors. Each carries its measured value as a load-bearing comment with a
`BENCHMARKS.md` pointer; a crossover changed without a re-measurement is a
regression with a green test suite.

## Why compose `butterfly-fft` rather than reimplement

Fast multipoint evaluation over a size-2^k additive subspace **is** the additive
FFT. `butterfly-fft` already owns that layout, its novel/Cantor bases, its
transform buffers, and the only SIMD `unsafe` in the L1 tier
(`butterfly-fft/src/core/transform.rs:307-680`; basis convert at
`basis/convert.rs:42-183`). Reimplementing it here would:

- **Duplicate the one hard, `unsafe`-carrying kernel** in a crate that
  `#![forbid(unsafe_code)]` (U4). The transform's vectorization is precisely the
  code that should live once, behind a sealed bound, and be audited once.
- **Fork the basis contract.** A monomial-basis `Polynomial` must convert to the
  novel basis before a `TransformPlan` and back after; owning a second copy of
  `monomial↔novel` would let the two drift, and the transform's correctness
  depends on that conversion matching the plan exactly.
- **Break the object rule.** The transform is `butterfly-fft`'s object; the ring
  is ours (ground rule 1). A crate that owns both owns two objects.

So U2 is a boundary, not a preference: `eval/transform.rs` calls
`TransformPlan::forward/inverse`, `ShiftedPlan` for cosets, and
`monomial_to_novel`/`novel_to_monomial`, and `univariate` reimplements only the
arbitrary-point path (`eval/horner.rs`, `eval/multipoint.rs`) that no transform
covers. The `fft` feature is default-on so the structured paths are present for
`gs-engine`/`hasse`/`funcfield`/`reed-muller`; a transform-free consumer
(`syndrome-engine`) builds `--no-default-features`, and the domain type stays
total by falling back to the arbitrary path over the enumerated subspace. This is
the mirror of the `fgf` composition rule (U1): the kernel owner is upstream, and
`univariate` is a consumer of it, never a second implementation.

# Charter

`univariate` is the univariate-polynomial-ring node of the FEC stack: arithmetic
in F[x] over GF(2^m), underneath every consumer that has to multiply, divide,
gcd, evaluate, interpolate, find roots, or solve a key equation. Everything in
this planning set defers to this document.

## What `univariate` is

`univariate` is a **ring, not a decoder.** It owns the dense monomial-basis
polynomial type over a binary field, every ring operation on it, division and
gcd with cofactors, the truncated/partial EEA that reconstructs a rational
function, root-finding (Chien, equal-degree factorization, power-series lifting),
evaluation and interpolation over arbitrary and structured point sets, and
truncated power-series inversion.

It exists because a full, `fgf`-dispatched, packed-byte univariate library
already lives in this repository — trapped as the private `poly/` and `roots/`
modules of `gs-engine` — and cannot be reused by the five other crates that need
it. `univariate` is the extraction; `gs-engine` becomes its first consumer. The
classical decoder primitives those consumers also need — cofactor EEA, truncated
EEA, Chien search, power-series inverse — exist nowhere in the tree, and this
crate is where they are built once.

## The one-sentence scope

**Given polynomials over GF(2^m), compute in the ring — multiply, divide, gcd,
evaluate, find roots, interpolate, and work modulo `x^t` — and never construct a
code, a matrix, or a transform buffer.**

## Scope boundary

Copy this into `src/lib.rs` and `AGENTS.md` verbatim:

> Given polynomials over GF(2^m), compute in the ring — multiply, divide, gcd,
> evaluate, find roots, interpolate, and work modulo x^t, and never construct a
> code, a matrix, or a transform buffer. Field arithmetic comes from `fgf`;
> transform-domain evaluation/interpolation from `butterfly-fft`; matrices from
> `gfm`. This crate receives polynomials and returns polynomials.

This is not decoration. It is the rule that decides every borderline item in
[`01-landscape.md`](01-landscape.md).

### In scope

| Area | Why it is here |
| --- | --- |
| **`Polynomial<F: FieldKernels>`** — dense monomial-basis, packed `Vec<u8>` LE (`F::BYTES`/coeff), empty = zero, `normalize()` trims high zeros | The single core type. It is the object the whole crate is about; every operation returns or consumes it. Lifted from `gs-engine/src/poly/mod.rs:62`. |
| Ring ops: `add`/`add_assign`/`add_scaled(_assign/_shifted)`, `scale(_assign)`, `shifted` (X^n·), `multiply_x_plus`, schoolbook `multiply`, `multiply_truncated`, char-2 `square_into` (O(deg)), `formal_derivative`, `hasse_derivative(order)` + `evaluate_hasse`, `binomial_odd` | The whole ring API. Squaring is a bit-spread in characteristic 2, not a multiply. Lifted from `gs-engine/src/poly/arithmetic.rs:12-263,434-571`. |
| Division: `div_rem(_into)`, `exact_divide`, `remainder`, `monic`, `multiply_mod`/`square_mod`/`pow_mod`, `divide_by_x_power(_into)`/`x_valuation` | The quotient-ring operations consumers build modular arithmetic on. Lifted from `arithmetic.rs:203-347`. |
| GCD & cofactors: existing plain monic `gcd`, PLUS the **new** extended Euclidean with Bézout cofactors, PLUS the **new** truncated/partial EEA (early stop at remainder degree < t) | The reusable key-equation / Padé primitive. Plain gcd exists (`arithmetic.rs:271-281`); the cofactor and truncated forms are a GAP built fresh (§01). |
| Root-finding: `base_field_roots` (gcd(p, X^\|F\|+X) + trace/Cantor–Zassenhaus splitting), `roth_ruckenstein` and `alekhnovich` power-series lifting, PLUS the **new** classical **Chien search** and a standalone linearized/affine root solver | Chien is the cheap O(n·deg) scan for the small locators of bounded-distance decoding; equal-degree wins for larger-degree factor extraction. Both ship. Extract `gs-engine/src/roots/field_roots.rs:104-330`, `roots/{roth_ruckenstein,alekhnovich}.rs`; Chien + linearized are a GAP. |
| Evaluation: Horner `evaluate`/`evaluate_many`; the **new** subproduct-tree multipoint evaluation over ARBITRARY point sets; structured subspace/coset evaluation by composing `butterfly-fft` | Horner and small point sets exist (`arithmetic.rs:138-160`); the arbitrary-point subproduct tree is a GAP; the structured path composes the transform (U2). |
| Interpolation: Newton basis + incremental Newton interpolation; Lagrange (arbitrary points, subproduct tree); inverse-FFT interpolation via `butterfly-fft` | Unifies the THREE duplicate Lagranges the ecosystem grew (§01). Lift `gs-engine/src/interpolation/plan.rs` `build_newton_basis`, `module.rs:470-484` `interpolate_newton_into`. |
| Truncated power series: **new** `inverse mod x^t` (Newton iteration), power-series division / reciprocal, `mul_trunc` | The interleaved-RS and Padé consumers need series inversion; `multiply_truncated` already exists, the inverse is a GAP built fresh. |
| FFT-domain product: `multiply_batch_truncated` via `butterfly-fft`, plus a **new** Karatsuba middle tier | The product hierarchy is schoolbook ↔ Karatsuba ↔ AFFT; extract the AFFT tier (`gs-engine/src/poly/afft.rs:144-353`), Karatsuba is a GAP (only schoolbook and AFFT exist today). |
| **`EvaluationDomain<F>`** (arbitrary / additive-subspace / affine-coset) + backend selection | The dispatch-policy home. Lifted from `gs-engine/src/domain.rs`, `cost.rs`. |

### Out of scope

| Area | Where it belongs |
| --- | --- |
| Field arithmetic, SIMD field kernels, backend detection | `fgf`. No field loops, no `unsafe` here. |
| The additive-FFT transform itself, novel/Cantor bases, transform buffers | `butterfly-fft`. `univariate` composes, never reimplements. |
| Matrices, Gaussian elimination, Cauchy/Vandermonde matrix inversion, polynomial-**matrix** row reduction / weak-Popov (`gfm/src/poly.rs`) | `gfm`. `univariate` owns polynomials, not matrices of polynomials. |
| **Bivariate** polynomials (`gs-engine/src/poly/bivariate.rs`) and the GS interpolation / list-decoding machinery | `gs-engine`. `univariate` is strictly univariate; the bivariate layer is a consumer of it. |
| Syndromes, key-equation *decoding orchestration*, Chien-over-received-word + Forney | `syndrome-engine`. `univariate` supplies the primitives (EEA, Chien, eval, derivative); it does not know what a received word is. |
| Wire formats, code parameters, prime fields GF(p) | Consumers; and `fgf` cannot express GF(p). |

### Deliberately deferred, not rejected

- **Half-GCD / fast EEA (Brent–Gustavson–Yun).** The constant-factor
  BM/schoolbook EEA wins at RS/BCH sizes (degrees in the hundreds); `gfm`'s own
  plan reaches the same conclusion (`gfm/.plans/03-algorithms.md:296-323`). Gate
  on a consumer at degree ≳ thousands.
- **Frobenius-based fast root-finding** beyond what `gs-engine` already ships.
- **Thread parallelism (`rayon`).** Feature-gated, default off, symbol/batch axis
  only.

## Load-bearing invariants

Numbered because the roadmap's acceptance criteria cite them. Violating one is a
bug even when the tests pass.

**U1 — Compose `fgf`, never re-host.** Every coefficient operation is
`fgf::field::Elem` (scalar) or `fgf::ops::*` (packed), dispatched above the
lane-bytes crossover. No hand-rolled field loop; the packed-vs-scalar crossover
is measured (`gs-engine/src/poly/arithmetic.rs:573`), not guessed.

**U2 — Compose `butterfly-fft` for structured domains.** Subspace/coset
multipoint eval and interpolation call `TransformPlan` + `monomial↔novel`; the
crate reimplements only the arbitrary-point (Horner / subproduct-tree) path.
There is no second additive-FFT.

**U3 — One polynomial type, packed and canonical.** `Polynomial<F>` stores
`Vec<u8>` packed LE, always normalized (no trailing-zero coefficients); the zero
polynomial is the empty buffer. The degree of zero is a documented sentinel.

**U4 — No `unsafe`.** `#![forbid(unsafe_code)]`. The only SIMD is upstream in
`fgf` / `butterfly-fft`.

**U5 — Steady-state zero allocation.** Every hot op has a `*_into` /
scratch-owning form; proven by `tests/zero_alloc.rs` under a counting global
allocator. `gs-engine` already follows this (`_into` + `*Scratch`).

**U6 — `inv(0) == 0` is inherited; test explicitly.** Division and root logic
never infer "not a root" or "invertible" from a zero result; leading
coefficients and pivots are checked with `is_zero()`.

**U7 — Measured crossovers are load-bearing comments.** Schoolbook ↔ Karatsuba ↔
AFFT product, Chien ↔ equal-degree root-finding, Newton ↔ transform
interpolation: each dispatch carries its measured threshold with a
`BENCHMARKS.md` pointer. Reuse `gs-engine`'s existing constants
(`AFFT_*_CROSSOVER`, `MODULE_INTERPOLATION_CROSSOVER=8`,
`DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER=20000`, `cost.rs` selectors) rather than
re-deriving.

**U8 — Determinism is a wire property where it feeds decoders.** Root
enumeration order and interpolation point order are documented and stable
(consumers map roots→positions). Changing them is a format break, not a refactor.

**U9 — The public boundary validates; kernels assert.** Geometry (lengths,
degrees, domain sizes) is checked with typed errors naming the offending value
*and* the limit at every public entry; internal kernels use `debug_assert!`.

## Settled decisions

1. **`univariate` takes the `poly` node; the graph is updated in both notations.**
   Ground rule 1 names "polynomial rings (`poly`)" as a sanctioned object — this
   crate IS that object, renamed. README mermaid + D2 and the AGENTS crate index
   change `poly` → `univariate` and every `poly-->X` edge → `univariate-->X`. See
   [`09-poly-rename.md`](09-poly-rename.md).
2. **The truncated-EEA / connection-polynomial primitive lives here, not `gfm`.**
   It is gcd-with-cofactors over F[x] — a polynomial-ring operation adjacent to
   `div_rem`/`gcd` this crate already owns. `gfm`'s *unshipped* Hankel/BM claim
   (`gfm/.plans/00-charter.md:49`, `02-architecture.md:50`,
   `03-algorithms.md:296`) is the linear-algebra / Wiedemann view over a scalar
   sequence; it must be re-scoped to that matrix use or compose `univariate`, so
   the primitive is written once. **This is an open cross-plan seam (see
   [`08-risks.md`](08-risks.md) D1); flagged, not silently overridden.** No code
   collision today — `gfm`'s `hankel.rs` does not exist in the tree (only
   `gfm/src/poly.rs` weak-Popov ships).
3. **`butterfly-fft` is an optional (default-on) `fft` feature dependency.** The
   core ring + gcd/EEA + division + Chien + Horner eval + power series need only
   `fgf`; a minimalist consumer (e.g. `syndrome-engine`, which needs no
   transform) builds `--no-default-features` to drop `butterfly-fft`. `gs-engine`
   / `hasse` / `funcfield` / `reed-muller` keep the default and get the FFT paths.
4. **Chien search AND equal-degree factorization both ship.** They win in
   different regimes: Chien is the cheap O(n·deg) scan for the small locators of
   bounded-distance decoding; the gcd(p, X^q+X)+split route wins for
   larger-degree factor extraction (`gs-engine`'s list-decoding roots). `cost.rs`
   picks. Do not delete `gs-engine`'s route in favour of Chien or vice versa.
5. **Bivariate stays in `gs-engine`.** `univariate` is univariate; extracting the
   bivariate layer would put a second object in the crate (ground rule 1).

## Dependencies

Stated in full in [`05-conventions.md`](05-conventions.md) § Dependencies. In one
line: **`fgf` required and pinned by rev, `butterfly-fft` optional under the
default-on `fft` feature, `rayon` optional and off, `criterion`+`proptest` for
tests/benches, nothing else.**

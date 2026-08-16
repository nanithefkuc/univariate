# Landscape: what already exists, and what this crate extracts

Working tree of 2026-08-16. Every `file:line` refers to that state.

Verdict legend:

- **TAKE** — the algorithm and roughly its shape move into `univariate`; the
  original becomes a caller or is deleted.
- **LIFT** — the univariate parts of a larger routine move; `univariate` owns the
  polynomial element ops and the original keeps the surrounding structure
  (matrix, module) and re-imports.
- **LEAVE** — stays where it is, with the reason recorded so it is not
  relitigated; a duplicate flagged for later migration is `LEAVE (flag)`.

Unlike `gfm`, which collapsed six independent implementations, this is a
**single-source extraction**: `gs-engine/src/poly/` + `gs-engine/src/roots/` is
one canonical, mature, `fgf`-dispatched library, and the rest of the ecosystem
grew thin duplicates (three Lagranges) around it because it was private.

## 1. The extraction inventory

| # | Item | Location (2026-08-16) | Shape | Verdict |
| - | ---- | -------------------- | ----- | ------- |
| 1 | `Polynomial<F>` type | `gs-engine/src/poly/mod.rs:56-243` | packed `Vec<u8>` LE, `normalize`, constructors/accessors, `PolynomialError{Config, DivisionByZero, NonExactDivision}` | **TAKE** — becomes the crate's core type; `gs-engine` re-imports. |
| 2 | Ring / field-poly arithmetic | `gs-engine/src/poly/arithmetic.rs:12-571` | add/scale/shift/mul (schoolbook)/mul_truncated/eval (Horner)/eval_many/hasse/formal_derivative/div_rem/exact_divide/monic/gcd/mul_mod/pow_mod/square_into/binomial_odd/`use_packed_kernel` crossover | **TAKE** — the whole ring API. |
| 3 | AFFT batched product | `gs-engine/src/poly/afft.rs:144-353` | `multiply_batch_truncated[_with]`, `PolynomialProductScratch`, `ProductStrategy{Auto, Schoolbook, Afft}`, crossovers | **TAKE (composes `butterfly-fft`)** — gated by the `fft` feature. |
| 4 | Base-field roots | `gs-engine/src/roots/field_roots.rs:104-330` | gcd(p, X^\|F\|+X)+trace-split (Cantor–Zassenhaus), reusable-scratch poly mul_mod/square_mod | **TAKE.** |
| 5 | Roth–Ruckenstein / Alekhnovich lifting | `gs-engine/src/roots/{roth_ruckenstein,alekhnovich}.rs` | power-series root lifting, `AlekhnovichScratch`, `DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER=20000` | **TAKE** — the power-series root backends. |
| 6 | Newton interpolation + vanishing-poly build | `gs-engine/src/interpolation/plan.rs`, `module.rs:470-484` | `build_newton_basis` (O(n²)), `interpolate_newton_into`, ∏(X−αᵢ) product/subspace | **LIFT** — the univariate parts move; the *module* (poly-matrix, weak-Popov) interpolation stays in `gs-engine` and keeps calling `gfm::weak_popov`. |
| 7 | `EvaluationDomain` + backend selectors | `gs-engine/src/domain.rs`, `cost.rs` | arbitrary/additive/coset domains, `EvaluationBackend`, product/root/interp selectors | **TAKE** — the dispatch-policy home. |
| 8 | Vandermonde-inverse hand-rolled Lagrange | `gfm/src/structured/vandermonde.rs:127-160` | master poly ∏(x+xₘ), synthetic division, Horner normalizer, on raw `Vec<F::Elem>` | **LEAVE (flag)** — `gfm` keeps the *matrix* type; its poly product/division/eval SHOULD re-base on `univariate` once it exists (consumer-migration phase). A duplicate to retire, not `univariate`'s code. |
| 9 | Exponent-domain Lagrange (systematic generator) | `srs/src/afft/generator.rs:11-45` | log/exp-domain Λ(x)/((x⊕d)·Λ'(d)) | **LEAVE** — `srs`'s transform-domain erasure path; a different representation, not extracted, but noted as a third Lagrange the ecosystem should not grow a fourth of. |
| 10 | Weak-Popov / poly-matrix row reduction | `gfm/src/poly.rs` (whole file) | `weak_popov*`, `WeakPopovRow/Basis`, `PopovLeadingTerm`, `WeakPopovScratch` | **LEAVE** — matrices of polynomials are `gfm`'s linear-algebra object; `univariate` supplies the F[x] *element* ops `gfm` composes, not the reduction. |
| 11 | Field-vector kernels | `fgf/src/ops.rs:322-774` | AXPY/scatter/gather/elementwise packed kernels | **LEAVE** — `fgf`'s, composed via U1. |
| 12 | `binomial_odd` duplicated inside `gs-engine` | `gs-engine/src/poly/arithmetic.rs:569` (+ historically interpolation) | Lucas/Sierpiński parity | **TAKE** — becomes one `univariate` helper; `gs-engine`'s internal duplicate collapses. |

**Count of univariate-polynomial libraries after this plan lands: one**
(`univariate`), plus `srs`'s deliberately different exponent-domain generator and
`gfm`'s matrix element ops re-based on it.

## 2. Gaps — primitives that exist NOWHERE, build fresh

Every item below is named as a need by a consumer plan and implemented nowhere in
the tree.

| Gap | Evidence |
| --- | --- |
| **Extended Euclidean with Bézout cofactors** over GF(2^m)[x] | Only plain monic `gcd` exists today (`gs-engine/src/poly/arithmetic.rs:271-281`); no routine returns the cofactors `s, t` with `s·a + t·b == g`. |
| **Truncated / partial EEA** (key-equation / Padé primitive, early stop deg < t) | The reusable connection-polynomial solver `syndrome-engine` and `contort`'s interleaved-RS path need. Named nowhere in shipping code. |
| **Classical Chien search** (root search by domain scan) | `gs-engine` uses gcd+split for roots; the cheap O(n·deg) locator scan of bounded-distance decoding is absent. |
| **Standalone linearized / affine (2-linearized) root solver** | No affine-root routine exists outside the equal-degree machinery. |
| **Truncated power series**: `inverse mod x^t` (Newton), power-series division, reciprocal | `multiply_truncated` exists; the inverse and division do not. |
| **Karatsuba multiplication** (middle tier) | Only schoolbook and AFFT exist; the constant-factor middle tier between them is unfilled. |

## 3. Consumer demand table

What each downstream crate actually needs, from its own plans and code.

| Consumer | Needs | Status |
| --- | --- | --- |
| `syndrome-engine` (L2, planned) | truncated-EEA / BM key-solver primitive, Chien, Horner eval (syndromes), `formal_derivative` + eval + batch-inv (Forney), `div_rem`/`gcd`. **Primary driver.** Needs NO transform → uses `--no-default-features`. | Planned crate; no code. Drives the GAP primitives. |
| `gs-engine` (L2, impl) | Its own `Polynomial`/roots/interp — the extraction source; becomes the reference consumer post-migration. | Implemented; the extraction source. Migration is Phase 6. |
| `contort` (L3, impl) | Interleaved-RS defers common-locator / key-equation + Chien + per-row Forney upstream (`contort/.plans/03-algorithms.md:245-285`); gets them via `syndrome-engine`, which gets primitives here. | Implemented; consumes indirectly. |
| `hasse` (L1, planned) | `evaluate_hasse` / `hasse_derivative` + subspace multipoint eval (edges `poly→hasse`, `butterfly→hasse`). | Planned crate. |
| `funcfield` (L1, planned) | Ring ops + gcd/division (edge `poly→funcfield`). | Planned crate. |
| `reed-muller` (L3, planned) | Eval/interpolation over additive subspaces. | Planned crate. |
| `gfm` (L1, impl) | Vandermonde inverse's hand-rolled poly ops → migrate to `univariate` (item 8), retiring a duplicate. | Implemented; migration is a consumer-side follow-up. |

Read against the charter, one entry deserves emphasis: **`gs-engine`'s hot path
is the extraction source, not a new customer.** The value of this crate to
`gs-engine` is not speed — the code is already `gs-engine`'s — it is that four
other crates stop reimplementing the same ring. The primary *new* demand comes
from `syndrome-engine`, and it is demand for the GAP primitives (cofactor EEA,
Chien, power-series inverse) that no permissive Rust crate provides.

## Provenance

Single-source extraction. Reconnaissance covered `fgf`, `butterfly-fft`,
`gs-engine`, `gfm`, `srs`, `contort`, and the planned consumer set at the
2026-08-16 working tree; every `file:line` above is scout-verified at that state.
The extraction inventory (§1) and the gap list (§2) are the primary evidence for
the crate's existence: one mature private library that five crates cannot reach,
plus four decoder primitives that exist nowhere.

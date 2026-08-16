# Roadmap

Each phase ends on an observable, checkable result. A compile-only scaffold is
not evidence for an algorithm; the oracles named in
[`03-algorithms.md`](03-algorithms.md) are, and the load-bearing invariants
(**U1–U9**) named in [`00-charter.md`](00-charter.md) are cited by number.

**Phase 6 is the proof phase.** Phases 0–5 extract the ring, build the GAP
primitives, and wire the product tiers in isolation. The moment gs-engine — the
crate this ring was lifted out of — deletes its private `poly/` and `roots/` and
imports `univariate`, and its existing test suite passes unchanged, the
extraction is demonstrated rather than argued: the reference consumer's behavior
is the certificate. Nothing after Phase 6 is scope for the first release.

Phases 0–6 are release 0.1. Phase 7 is optional scope, gated on a consumer.

---

## Phase 0 — Scaffold and conventions

Manifest, lints, CI, `src/error.rs`, the AI-authorship `README.md`, `LICENSE`
(MIT), `CHANGELOG.md`, `AGENTS.md`, and the doc set per
[`05-conventions.md`](05-conventions.md). Delete the cargo template `add`
function.

`fgf` enters as a required dependency here, pinned by rev, because the layout
policy in Phase 1 needs `F::BYTES` and the packed-kernel dispatch. `butterfly-fft`
enters gated by the default-on `fft` feature; `rayon` enters gated by the
default-off `parallel` feature and is a no-op until Phase 7. Define
`PolynomialError` and `DomainError` with their variants (each carrying the
offending value *and* the limit, U9), and the pass-through `backend_for` if one
is needed. Nothing else.

**Acceptance:**

- `cargo fmt --all --check`; `cargo clippy --all-targets --all-features -- -D
  warnings` and `cargo clippy --all-targets --no-default-features -- -D
  warnings`;
- `cargo test`, `--all-features`, and `--no-default-features` (the last proving
  no `butterfly-fft` in the tree — R3);
- `RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps`;
- cross-builds to `aarch64-unknown-linux-gnu` and `wasm32-unknown-unknown`,
  `--no-default-features`;
- `cargo +1.89.0 build --all-features`;
- `backend_for::<Gf8B>()` returns the same value as `fgf`'s does on the same
  host, and returns `Scalar` under `SIMD_BACKEND=scalar`.

---

## Phase 1 — Core ring (extraction)

`poly/dense.rs`, `poly/ring.rs`, `poly/divide.rs`, `poly/gcd.rs` (plain gcd
only), `eval/horner.rs`, `cost.rs` (product/eval selectors only).

`Polynomial<F: FieldKernels>` — packed `Vec<u8>` LE, `F::BYTES` per coefficient,
empty buffer = zero, `normalize()` trims trailing zeros, degree-of-zero a
documented sentinel (U3). Lifted from `gs-engine/src/poly/mod.rs:56-243`. Ring
ops (`add`/`add_assign`/`add_scaled(_assign/_shifted)`, `scale(_assign)`,
`shifted`, `multiply_x_plus`, schoolbook `multiply`, `multiply_truncated`,
char-2 `square_into`, `formal_derivative`, `hasse_derivative`/`evaluate_hasse`,
`binomial_odd`) from `arithmetic.rs:12-263,434-571`. Division (`div_rem(_into)`,
`exact_divide`, `remainder`, `monic`, `*_mod`, `divide_by_x_power`/`x_valuation`)
from `arithmetic.rs:203-347`. Plain monic gcd from `arithmetic.rs:271-281`.
Horner `evaluate`/`evaluate_many` from `arithmetic.rs:138-160`. Every hot op
lands with its `_into` form (U5). The packed-vs-scalar crossover comment is
reused, not re-derived (`arithmetic.rs:573`, U7). No hand-rolled field loop; all
coefficient work is `fgf::field::Elem` or `fgf::ops::*` (U1).

**Acceptance:**

- §03 oracles 1 (product tiers, schoolbook only this phase — byte-identical
  against the naive convolution), 2 (`div_rem`: `q·b + r == a`, `deg r < deg b`),
  and 3 (plain gcd: `g | a`, `g | b`, `gcd(a,0)`/`gcd(a,a)` cases) pass on
  randomized inputs at every degree including the empty/zero and unit-divisor
  cases;
- §03 oracle 5 (Horner) passes: `evaluate` and `evaluate_many` against a naive
  reference, and the char-2 `square_into` agrees with `multiply(p, p)`;
- all pass across `Gf8B`, `Gf8D`, `Gf16`, `Gf32`, `Gf64`, and `FanPaar8/16/32/64`;
- `Polynomial` is always canonical: no trailing-zero coefficient survives any
  constructor or mutating op, checked as a property (U3);
- `inv(0)==0` is exercised: a zero leading coefficient / zero divisor is caught
  by `is_zero()` and returns `DivisionByZero`, never a silent wrong quotient
  (U6);
- `tests/zero_alloc.rs` passes under the counting global allocator for every
  `*_into` path shipped this phase (U5).

---

## Phase 2 — Extended gcd, truncated EEA, power series (GAPS)

`poly/gcd.rs` (extended + truncated EEA), `poly/series.rs`.

Extended gcd with Bézout cofactors (`s·a + t·b == g`), the truncated / partial
EEA (early stop at remainder degree `< t`, the key-equation / Padé primitive),
and truncated power series: `inverse_mod_x_t` by Newton doubling, series
division, reciprocal. These exist nowhere in the tree and are built fresh. The
truncated EEA is written *once* here; the gfm Hankel/BM seam is reconciled before
gfm implements its own (settled decision #2, [`08-risks.md`](08-risks.md) D1).

**Acceptance:**

- §03 oracle 3 extended: Bézout identity `s·a + t·b == g` by the independent
  multiply/add, `g | a` and `g | b` by `div_rem`, and the cofactor degree bounds
  asserted exactly;
- §03 oracle 4: the Padé identity (`t/s ≡ S mod x^t` by truncated multiply) and
  the **Dornstetter cross-check** against a reference Berlekamp–Massey written
  once in `tests/` — same connection polynomial on the same sequence, plus the
  LFSR-reproduction property;
- §03 oracle 9: `a · a_inv ≡ 1 (mod x^t)` by truncated multiply, and Newton
  doubling agrees with the linear schoolbook series inverse; a zero constant
  term returns the documented error, not a silent zero (U6);
- cross-checks against `galois` (`gcd`, `berlekamp_massey`) on fixed-seed
  vectors;
- `tests/zero_alloc.rs` covers the new `*_into` forms (U5).

---

## Phase 3 — Root-finding

`roots/equal_degree.rs`, `roots/lift.rs` (extract), `roots/chien.rs`,
`roots/linearized.rs` (GAP), `cost.rs` (root selector).

Extract `base_field_roots` (gcd(p, X^{|F|}+X) + Cantor–Zassenhaus split) from
`roots/field_roots.rs:104-330` and the Roth–Ruckenstein / Alekhnovich lifting
backends (with `AlekhnovichScratch`, `DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER=20000`)
from `roots/{roth_ruckenstein,alekhnovich}.rs`. Add the classical Chien search
and the standalone linearized / affine (2-linearized) root solver. `cost.rs`
selects Chien vs equal-degree; both ship, neither replaces the other (settled
decision #4). Enumeration order is frozen (U8).

**Acceptance:**

- §03 oracles 7 and 8: every returned root satisfies `evaluate(p, r) == 0` by
  the independent Horner; Chien and equal-degree agree on root **sets** and
  **counts**, and both respect the degree bound;
- the two lift backends (Roth–Ruckenstein, Alekhnovich) return the same root set
  on the same input;
- the linearized solver's roots match the Chien enumeration on affine inputs;
- the `cost.rs` selector picks the cheaper backend on measured sizes, with the
  crossover recorded in `BENCHMARKS.md` (U7);
- root enumeration order is stable across runs and backends — a frozen-order
  test (U8, R4);
- cross-checks against `galois` `.roots()` on fixed vectors;
- `tests/zero_alloc.rs` covers the scratch-owning root paths (U5).

---

## Phase 4 — Evaluation domains and interpolation

`eval/multipoint.rs`, `eval/newton.rs`, `eval/domain.rs`, `eval/transform.rs`
(`#[cfg(feature="fft")]`).

Extract Newton-basis interpolation (`build_newton_basis`,
`interpolate_newton_into`, `MODULE_INTERPOLATION_CROSSOVER=8`) from
`interpolation/plan.rs` and `module.rs:470-484` — the univariate parts only, the
module/weak-Popov interpolation stays in gs-engine (item 6, LIFT). Add the
arbitrary-point subproduct-tree multipoint evaluation and Lagrange
interpolation. `EvaluationDomain<F>` (arbitrary / additive-subspace /
affine-coset) with backend selection from `domain.rs`/`cost.rs`. Compose
`butterfly-fft` under the `fft` feature for subspace `forward`/`inverse` +
`monomial↔novel` (U2).

**Acceptance:**

- §03 oracle 5 (multipoint): subproduct-tree evaluation equals per-point Horner
  on arbitrary point sets, and equals `butterfly-fft` `forward` over a subspace;
- §03 oracle 6 (interpolation): round-trip (interpolate then `evaluate_many`
  returns the values, `deg < n`); Newton, Lagrange, and inverse-FFT **agree** on
  their shared domains — the differential the ecosystem's three duplicate
  Lagranges never had (gfm Vandermonde and srs generator cited as the duplicates
  to reconcile, items 8 and 9);
- the `fft`-gated paths are absent from the `--no-default-features` API surface,
  asserted by a build without the feature (R3);
- domain-geometry errors (`DomainError::LengthMismatch`, `NotSubspace`) name the
  value and the limit (U9);
- cross-checks against `galois` interpolation on fixed vectors;
- `tests/zero_alloc.rs` covers the `_into` interpolation/eval paths (U5).

---

## Phase 5 — Product tiers and Karatsuba

`poly/karatsuba.rs`, `poly/afft.rs` (extract), `cost.rs` (product crossovers).

Extract the AFFT batched product (`multiply_batch_truncated[_with]`,
`PolynomialProductScratch`, `ProductStrategy{Auto,Schoolbook,Afft}`,
`AFFT_*_CROSSOVER`) from `poly/afft.rs:144-353`, gated by `fft`. Add the
Karatsuba middle tier (GAP). Wire the measured schoolbook↔Karatsuba↔AFFT
crossovers through `cost.rs` (U7).

**Acceptance:**

- §03 oracle 1 in full: the three tiers are **byte-identical** on every operand
  pair over their overlapping degrees — Karatsuba vs schoolbook everywhere, AFFT
  vs both where the subspace covers the product degree;
- the crossovers are measured, not guessed, and recorded in `BENCHMARKS.md` with
  hardware, rustc version, and command line following `fgf`'s measurement
  hygiene (U7);
- the AFFT path is `fft`-gated and absent under `--no-default-features` (R3);
- no regression on any Phase 1–4 benchmark;
- `tests/zero_alloc.rs` covers `PolynomialProductScratch` reuse (U5).

---

## Phase 6 — Consumer migration (clean cutover). **This decides whether 0.1 exists.**

The extraction's proof. `gs-engine` deletes its private `poly/` and `roots/` and
imports `univariate`; the univariate parts of `gs-engine/src/interpolation/` are
deleted, the module/weak-Popov reduction kept (item 6). `gfm`'s Vandermonde
inverse re-bases its hand-rolled poly product/division/eval on `univariate`
(item 8), retiring that duplicate. No `#[deprecated]`, no re-export shim, no
`compat` module anywhere.

**Acceptance — this decides whether the extraction is correct:**

- **gs-engine's full test suite passes unchanged** after substitution — the
  tests were written against behavior, so a passing suite is the substitution's
  proof, and gs-engine's existing tests *are* the differential oracle (R1, no
  dual-maintenance window);
- gfm's suite passes unchanged after its Vandermonde poly ops re-base on
  `univariate`;
- root enumeration and interpolation-point order are byte-stable across the
  cutover — the frozen-order property consumers' position maps depend on (U8,
  R4);
- each migrated crate's benchmarks are re-run and the before/after ratio
  recorded; a regression greater than 5% blocks that crate's migration and is
  investigated, not accepted;
- a dependency-tree check confirms `univariate` at default features carries
  `fgf` and `butterfly-fft`, and that a `--no-default-features` build carries
  `fgf` alone — the seam `syndrome-engine` relies on (R3).

---

## Phase 7 — Parallelism (optional)

`parallel` feature, `rayon`, default off. Symbol/batch axis only: `par_iter` over
independent columns of `multiply_batch_truncated`, over independent point-set
leaves of multipoint evaluation, and over independent rows of a batched
interpolation — never inside a single polynomial's coefficient loop, which is
inherently sequential.

None of this starts without a consumer that demonstrably needs it.

**Acceptance:**

- identical results with the feature on and off, on every fixture;
- scaling measured on 2/4/8/16 threads and recorded in `BENCHMARKS.md`;
- the crate still builds and tests `--no-default-features` and on
  `wasm32-wasip1`;
- no `rayon` in the dependency tree without the `parallel` feature, asserted by
  a CI dependency-tree check.

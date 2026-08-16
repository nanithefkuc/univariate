# Conventions

Distilled from `fgf/AGENTS.md`, `gfm/.plans/05-conventions.md`, and
`gs-engine/src/lib.rs`. Where the siblings disagree, the choice and its reason
are stated.

## Dependencies

**Two runtime dependencies, one optional. This is a hard constraint, not a
preference.**

| Dependency | Kind | Version / pin | Why | Feature-gated |
| --- | --- | --- | --- | --- |
| `fgf` | **runtime, required** | git, **pinned by rev** | Field arithmetic (`field::Elem`/`Field`), the packed-buffer kernels (`ops::*`), and backend dispatch (`kernel::{Backend, FieldKernels, backend_for}`). The U1/U2 seam forbids re-hosting any of it. | no |
| `butterfly-fft` | runtime, optional | git, **pinned by rev** | Additive-FFT multipoint evaluation/interpolation over additive subspaces and affine cosets, plus monomial↔novel basis conversion. `univariate` composes it (U2); it never reimplements the transform. | `fft`, **default on** |
| `rayon` | runtime, optional | `1` | Batch/symbol-axis parallelism, Phase 7. | `parallel`, **default off** |
| `criterion` | dev | `0.8`, `default-features = false`, `features = ["cargo_bench_support"]` | Benchmarks. Matches `gfm`. | — |
| `proptest` | dev | `1`, default features | Property tests, fixed seed. | — |

**Nothing else. Ever.** Specifically ruled out:

- **`gfm`** — matrices. The polynomial ring does not eliminate, factor matrices,
  or solve a linear system; it has no need for `Ple`, Cauchy/Vandermonde inversion,
  or weak-Popov reduction. The polynomial-*matrix* row reduction and the
  matrix-form Cauchy/Vandermonde inverse are `gfm`'s object (ground rule 1);
  `univariate` supplies the F[x] *element* ops `gfm` composes, not the reduction.
  An edge here would put a second object in the crate.
- **`sgraph`** — sparse graphs and peeling; a different object entirely.
- **`lattica`** — exact integer arithmetic in a different ring with an
  overflow-based failure model that has no meaning over a finite field. We
  parallel its naming where operations correspond; we do not depend on it.
- **`syndrome-engine`** — the consumer; an edge would invert the layering
  (`syndrome-engine --> univariate`, not the reverse).
- **`archmage` / `simdispatch`** — we consume the re-export through `fgf`
  (`fgf::Backend`, `fgf::backend_for`) and never re-probe target features. No
  per-crate env override; the one stack-wide downgrade-only variable is
  `SIMD_BACKEND`, owned by `simdispatch`.
- **`thiserror`, `anyhow`** — hand-rolled errors. `fgf` and `gfm` both do this; a
  fourth error convention in the stack is worse than the boilerplate.
- **`nalgebra`, `ndarray`** — structurally closed (bounded `T: ComplexField`,
  needs norms/ordering); a GF(2^m) element cannot implement them.
- **`num-traits`, `zkcrypto/ff`, `ark-ff`** — prime-field trait ecosystems;
  `ff`'s `#[derive(PrimeField)]` cannot express GF(2^m).
- **A random-number crate** — tests use the fixed-seed LCG convention `fgf`
  already established (`noise(len, seed)`). No nondeterministic randomness, in
  tests or benches.

### Why `fgf` is pinned by rev

`gfm/Cargo.toml:34` pins by rev with the comment "a floating dependency is a
format-break risk", and it is right here for a stronger reason: `fgf`'s
`Backend` declaration order encodes capability and `Backend` is
`#[non_exhaustive]`; reordering variants is, in `fgf`'s own words, "a behavioral
and safety change". A floating git dependency turns any of that into a silent
break at the next unrelated rebuild. Pin by rev; swap the rev with every
ecosystem consumer in the same sitting.

Note: `gfm/.plans/09-fgf-rename.md` records that the field crate is `fgf` —
the pre-rename name is dead and gone from the tree. Every citation here uses
`fgf`; treat any stale pre-rename path in a sibling plan as `fgf::`.

### What we deliberately do *not* take from `fgf`

The **`internals` feature.** `gs-engine`'s `Cargo.toml` sets the precedent: it
declares `fgf` with no features and consumes only the stable surface. Under
`internals`, `fgf` exposes `kernel::xor`, internal matrix traits, and slab
helpers — none load-bearing for a polynomial crate. Taking it would couple us to
an explicitly non-semver surface for no capability we cannot otherwise reach.

`univariate` declares its **own** `internals` feature, gating its own unstable
types for benchmarks and downstream experimentation — exactly as `gfm` does.

## Cargo.toml

```toml
[package]
name = "univariate"
version = "0.0.0"
edition = "2024"
rust-version = "1.89"
license = "MIT"
publish = false
description = "Univariate polynomial arithmetic over GF(2^m): ring ops, gcd/EEA, division, root-finding, evaluation, interpolation, truncated power series."
categories = ["algorithms", "mathematics", "no-std"]
exclude = ["/.github", "/.plans"]

[features]
default  = ["std", "simd", "fft"]
std      = ["fgf/std"]
simd     = ["std", "fgf/simd"]          # simd implies std, as in fgf and gfm
fft      = ["std", "butterfly-fft/std", "butterfly-fft/simd"]
parallel = ["std", "dep:rayon"]
internals = []                          # our own unstable surface, not fgf's

[dependencies]
fgf            = { git = "https://github.com/nanithefkuc/fgf", rev = "<pin at Phase 0>" }
butterfly-fft = { git = "https://github.com/nanithefkuc/butterfly-fft", rev = "<pin at Phase 0>", optional = true }
rayon          = { version = "1", optional = true }

[dev-dependencies]
criterion = { version = "0.8", default-features = false, features = ["cargo_bench_support"] }
proptest  = { version = "1", default-features = true }

[[bench]]
name = "poly"
harness = false

[profile.bench]
lto = "thin"
codegen-units = 1

[package.metadata.docs.rs]
all-features = true
rustdoc-args = ["--cfg", "docsrs"]
```

Notes on the choices:

- **`publish = false`.** Nothing in this stack publishes until `fgf` reaches
  1.0.0. Reserve the name when the rest of the stack does, not before.
- **`edition = 2024`, `rust-version = "1.89"`.** Matches `fgf`, `gfm`, and
  `gs-engine` exactly. CI selects `+1.89.0` for the MSRV job.
- **`simd` implies `std`.** `fgf`'s backend cache is a `LazyLock`; a
  `simd`-without-`std` combination would compile and silently run scalar.
- **`fft` implies `std` and is in `default`.** The structured-domain
  eval/interp path composes `butterfly-fft`, which carries its own SIMD `unsafe`.
  Most L1 consumers want it; `syndrome-engine` builds
  `--no-default-features` to drop it. `fft` is the one feature a hot decoder
  opts out of.
- **`parallel` implies `std` and is not in `default`.** The crate stays usable
  on `wasm32-wasip1`.

## `src/lib.rs` header

```rust
//! Univariate polynomial arithmetic over GF(2^m).
//!
//! <scope boundary blockquote from 00-charter.md, verbatim>

#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]
#![warn(missing_docs, missing_debug_implementations)]
#![warn(clippy::pedantic)]
#![allow(clippy::cast_possible_truncation)] // <justification comment, required>
#![allow(clippy::module_name_repetitions)]

extern crate alloc;
```

`forbid`, not `deny` — following `gfm` and `gs-engine`: this crate writes no
intrinsics and composes no `unsafe` surface, so there is one implementation to
audit and it is upstream in `fgf`/`butterfly-fft`. (`butterfly-fft` is the one
L1 crate that carries SIMD `unsafe` and therefore uses `warn`; a *consumer* of
it still forbids.)

## Naming

Parallel to `gfm` and `gs-engine` wherever the operation corresponds, so a reader
who knows one knows the other:

| `gs-engine` (current private) | `univariate` |
| --- | --- |
| `Polynomial<F>` | same |
| `from_coefficients` / `constant` / `one` / `zero` | same |
| `add` / `add_scaled` / `scale` / `shifted` | same |
| `multiply` / `multiply_truncated` / `multiply_mod` | same |
| `div_rem` / `exact_divide` / `remainder` / `monic` | same |
| `gcd` | same; **`gcd_ext`** (new, cofactors) / **`truncated_eea`** (new, key-equation primitive) |
| `evaluate` / `evaluate_many` / `evaluate_hasse` | same |
| `hasse_derivative` / `formal_derivative` | same |
| `chien` (new) / `base_field_roots` / `roth_ruckenstein` / `alekhnovich` | same |
| `EvaluationDomain` / `EvaluationBackend` | same |

From `fgf/AGENTS.md`: conversions are `from_raw`/`to_raw`, `from_bytes`/`to_bytes`;
algebraic constants are uppercase; small hot wrappers are `#[inline]`; query APIs
returning a maybe-absent thing return `Option`.

From `srs` (carried into `gs-engine`): anything that returns a collection on a
per-symbol path also has an `_into(&mut …)` form, and the scratch object is
named `<Op>Scratch` (`PolynomialProductScratch`, `AlekhnovichScratch`,
`DecodeScratch`).

## Module and doc style

`lib.rs` and every `mod.rs` hold declarations only — no function bodies, no
`impl` blocks. All public items are re-exported at the crate root. Test
placement follows visibility: in-module `#[cfg(test)] mod tests;` for anything
touching private state, `tests/` for the public surface.

Hard bans in rustdoc: development history, phase numbers, downstream-private
crate names, and **throughput numbers**. `BENCHMARKS.md` is the only place that
holds a number, a CPU name, a rustc version, or a command line — the policy `fgf`
adopted. A doc comment states the decision and the mechanism, then points at
`BENCHMARKS.md`.

## Errors

Hand-rolled, one `src/error.rs`, small enums per failure domain,
`#[derive(Debug, Clone, PartialEq, Eq)] #[non_exhaustive]`, struct variants whose
fields carry both the offending value and the limit, manual `Display` with
inline-captured args, `std::error::Error` under `std`, and a `/// # Errors`
section on every fallible public function.

Two enums:

```rust
pub enum PolynomialError {
    Config { /* … */ },
    DivisionByZero { dividend_degree: Option<usize> },
    NonExactDivision { remainder_degree: usize },
    Degree { value: usize, limit: usize },
}
pub enum DomainError {
    LengthMismatch { expected: usize, got: usize },
    NotSubspace { size: usize },
}
```

## Tests

Four layers.

1. **`tests/oracles.rs`** — the naive reference implementations, written once
   and never optimized: schoolbook convolution, coefficient-by-coefficient
   Horner, long-division-by-hand, a plain extended-Euclid, a full-domain Chien
   scan, a linear Newton-interpolation loop. These are the oracles the suite leans
   on and they are deliberately slow.
2. **`tests/poly.rs`, `tests/roots.rs`, `tests/eval.rs`, `tests/series.rs`** —
   the public surface against those oracles, at the shapes and degrees the
   roadmap's acceptance criteria enumerate, across every field `fgf` exposes.
3. **`tests/zero_alloc.rs`** — an isolated integration binary with a counting
   global allocator, proving U5 for every `*_into` / scratch-owning path.
4. **`tests/foreign.rs`** — differential against **NTL** `GF2EX` (poly
   arithmetic, `MinPolySeq` for the EEA/BM equivalence, factoring for
   root-finding), gated on the library being present. **Skips loudly** — prints
   what it could not find and returns, never `#[ignore]`d, following `fgf`'s
   convention for hardware-gated tests.

Discipline:

- **Never use the implementation under test as its own oracle.** Where the
  natural oracle is a slower version of the same algorithm, a second,
  structurally different check is required — which is why the AFFT product path
  is checked against schoolbook *and* Karatsuba, and the EEA path against BM
  (Dornstetter) rather than against a slow EEA.
- A bug fix ships a failing regression first.
- Exact values, not predicates. `assert_eq!(deg(remainder), 3)`, never
  `assert!(remainder.degree() < 10)`.
- Fixed-seed LCG only; reuse `fgf`'s `noise(len, seed)` shape. No `rand`.
- Boundary shapes deliberately straddle every lane and unroll boundary.
- The zero polynomial, degree-0 constants, and a single point are cases, not
  edge cases.
- Panic tests use `#[should_panic(expected = "stable message fragment")]`.
- Every error path is tested for *state preservation*: snapshot before, provoke
  the error, snapshot after, assert equal.

## Benches

`benches/poly.rs`, criterion, `harness = false`. Groups: `multiply`
(schoolbook/Karatsuba/AFFT), `divide`, `gcd` (plain + extended + truncated EEA),
`eval` (Horner / subproduct-tree / butterfly), `roots` (Chien / equal-degree /
lifting), `series`.

Every case warms its scratch first — the numbers describe the allocation-free
steady state that `tests/zero_alloc.rs` proves. Baselines are deliberately not
committed.

Policy, from `gfm` and `fgf` jointly:

```sh
cargo bench --features internals -- --save-baseline before
# change
cargo bench --features internals -- --baseline before
```

**Do not land a performance change on the strength of reasoning alone.**
Measurement hygiene inherited from `fgf/BENCHMARKS.md`: interleave
`base, new, base, new, …`; take the **maximum** of at least three runs per key;
always keep an unchanged 1.00x control; treat 16–128 byte coefficient buffers as
noise, because dispatch-identical code varied by up to 1.36x there. Benchmarks
are not CI correctness checks.

## CI

`env: RUSTFLAGS: -D warnings`. Jobs, and the phase that adds each:

- **`check`** (Phase 0): `fmt --check`; `clippy --all-targets`; `clippy
  --all-targets --all-features`; `cargo test`, `--all-features`,
  `--no-default-features`; `cargo doc --all-features --no-deps` with
  `RUSTDOCFLAGS: -D warnings`. `Swatinem/rust-cache@v2`.
- **`cross`** (Phase 0): build-only for `aarch64-unknown-linux-gnu` and
  `wasm32-unknown-unknown`, `--no-default-features`.
- **`msrv`** (Phase 0): `dtolnay/rust-toolchain@1.89.0`, `cargo build
  --all-features`.
- **`backends`** (Phase 1): the full suite under `SIMD_BACKEND=v3_gfni_crypto`,
  `v3`, `v2`, and `scalar`, each asserting the backend it actually selected and
  failing loudly if the host could not provide it.
- **`no-fft`** (Phase 0): build and test `--no-default-features` (the
  `syndrome-engine` configuration); the `fft` paths must not leak into the core
  ring API.
- **`deps`** (Phase 0): a dependency-tree assertion. `cargo tree --depth 1` at
  default features must contain exactly `fgf` and `butterfly-fft`; at
  `--no-default-features` exactly `fgf`; and must never contain `gfm`,
  `sgraph`, `lattica`, `syndrome-engine`, or `archmage`.
- **`foreign`** (Phase 2, optional): installs NTL and runs `tests/foreign.rs`.
  Allowed to fail the *install*, not the comparison.

## Draft `AGENTS.md`

````markdown
# univariate

> Given polynomials over GF(2^m), compute in the ring — multiply, divide, gcd,
> evaluate, find roots, interpolate, and work modulo x^t, and never construct a
> code, a matrix, or a transform buffer. Field arithmetic comes from `fgf`;
> transform-domain evaluation/interpolation from `butterfly-fft`; matrices from
> `gfm`. This crate receives polynomials and returns polynomials.

## Non-negotiables

1. **Compose `fgf`, never re-host.** Every coefficient operation is
   `fgf::field::Elem` (scalar) or `fgf::ops::*` (packed), dispatched above the
   lane-bytes crossover. No hand-rolled field loop.
2. **Compose `butterfly-fft` for structured domains.** Subspace/coset
   multipoint eval and interpolation call `TransformPlan` + `monomial↔novel`; the
   crate reimplements only the arbitrary-point (Horner / subproduct-tree) path.
   No second additive-FFT.
3. **One polynomial type, packed and canonical.** `Polynomial<F>` stores `Vec<u8>`
   packed LE, always normalized; the zero polynomial is the empty buffer.
4. **No `unsafe`.** Forbidden at the crate root. The only SIMD is upstream.
5. **Steady-state zero allocation.** Every hot op has a `*_into` / scratch form;
   proven by `tests/zero_alloc.rs`.
6. **`inv(0) == 0` is inherited.** Division and root logic test pivots/roots with
   `is_zero()`; never infer "not a root" from a division result.
7. **Numbers live in `BENCHMARKS.md`.** Doc comments state the decision and the
   mechanism and point there.
8. **Do not land a performance change on reasoning alone.** A/B it, keep both
   twins compiled, record the ratio.
9. **Oracles stay independent.** An implementation is never its own test.

## Working here

- Edition 2024, MSRV 1.89. No toolchain pin; select `+1.89.0` for the MSRV job.
- Features: `default = ["std", "simd", "fft"]`; `simd` implies `std`; `fft`
  implies `std` and composes `butterfly-fft`; `parallel` is an off-by-default
  no-op placeholder; `internals` exposes this crate's unstable surface (never
  `fgf`'s — we do not enable it).
- `src/lib.rs` and every `mod.rs` hold declarations only — no function bodies,
  no `impl` blocks. Public items are re-exported at the crate root.
- Errors are hand-rolled in `src/error.rs`: small enums per failure domain,
  struct variants carrying the offending value and the limit, manual `Display`,
  `std::error::Error` under `std`. Every fallible public function documents
  `# Errors`.
- Test placement follows visibility: in-module `#[cfg(test)]` for private state,
  `tests/` for the public surface. Fixed-seed LCG only (`fgf`'s `noise(len, seed)`
  shape); no `rand`. Exact values, not predicates.
- The full check set:

  ```sh
  cargo fmt --all -- --check
  cargo clippy --all-targets --all-features -- -D warnings
  cargo clippy --all-targets --no-default-features -- -D warnings
  cargo test
  cargo test --all-features
  cargo test --no-default-features
  RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
  cargo build --target aarch64-unknown-linux-gnu --no-default-features
  cargo build --target wasm32-unknown-unknown --no-default-features
  cargo +1.89.0 build --all-features
  ```

- Benchmarks go through `criterion`; baselines are deliberately not committed.
  Measurement hygiene: interleave base/new, take the maximum of at least three
  runs, keep an unchanged 1.00x control, treat 16–128 byte buffers as noise.
- Commit subjects are at most ~10 words, shaped `univariate: short verb phrase`.
  What changed and why lives in the pull request and `CHANGELOG.md`.
````

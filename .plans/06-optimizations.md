# Optimizations

Ranked at the end by expected payoff ÷ implementation cost. Every entry names
the evidence. Where the evidence is a measurement `gs-engine` already made on
this ring, we reuse its threshold rather than re-deriving it (U7) — and where a
fast asymptotic algorithm was already weighed and *lost* at RS/BCH sizes, it is
recorded as deferred so it is not re-proposed.

Evidence tags: **[repo]** = measured or wired in this repository, cited by
`file:line`; **[pub]** = published result with a source in the
[bibliography](07-baselines.md); **[inf]** = inference from the above.

## The cost model, stated once

`univariate` owns **no field kernels** (U1). Every coefficient operation is a
`fgf::field::Elem` scalar or a `fgf::ops::*` packed-buffer call, dispatched above
a **measured** lane-bytes crossover, never a hand-rolled field loop. The field
layer below us is already tuned to the memory-bandwidth ceiling upstream in
`fgf`; a polynomial crate cannot make a GF(2^m) multiply faster and must not try.

So the wins here are never in a faster multiply. They are in three things, and
everything below is one of them or a rejection:

1. **Deleting coefficient operations** — asymptotically better algorithms that
   do fewer field ops for the same result: Newton-iteration series inversion
   `O(M(t))` against `O(t²)` schoolbook, subproduct-tree multipoint work
   `O(M(n) log n)` against `O(n²)`, and characteristic-2 squaring as an
   `O(deg)` bit-spread instead of a convolution.
2. **Reusing already-measured crossovers** (U7) — schoolbook↔Karatsuba↔AFFT,
   Chien↔equal-degree, Newton↔transform. The packed-vs-scalar crossover
   (`gs-engine/src/poly/arithmetic.rs:573` `use_packed_kernel`) **[repo]**, the
   AFFT crossovers (`poly/afft.rs`), `MODULE_INTERPOLATION_CROSSOVER = 8`, and
   `DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER = 20000` are load-bearing constants that
   already exist; re-deriving them is waste.
3. **Composing the transform, never re-hosting it** (U2) — the structured
   multipoint path *is* the additive FFT and belongs to `butterfly-fft`; the
   only speed we own on that path is the decision of *when* to enter it.

The corollary, stated so a future contributor does not re-open it: an
optimization that only makes the field multiply faster is out of scope by U1 —
it belongs upstream in `fgf`. An optimization that changes *how many* field
multiplies a polynomial operation issues is ours.

---

## 1. Packed-kernel dispatch above the lane-bytes crossover

**What.** Every coefficient-vector operation — `add_assign`, `add_scaled`
(AXPY), `scale`, elementwise product for evaluation lanes — routes to
`fgf::ops` packed-buffer kernels once the operand is wider than a measured
threshold, and stays a scalar `fgf::field::Elem` loop below it. This is U1 made
concrete: no polynomial-level field loop exists.

**Evidence.** The crossover is already measured and wired:
`gs-engine/src/poly/arithmetic.rs:573` `use_packed_kernel` decides scalar vs
packed on operand size **[repo]**. `fgf` exposes exactly the kernels the ring
needs — `add_assign::<F>(dst,src)` (`dst ^= src`), `mul_add::<F>` (AXPY) and its
`_with(PreparedCoefficient)` form, `mul_into`, `mul_assign`, `mul_elementwise`
for the evaluation/FFT lanes (`fgf/src/ops.rs:322-774`) **[repo]**. The buffers
are `F::BYTES`-packed `&[u8]`, which is the crate's canonical coefficient storage
(U3), so the dispatch is a length comparison, not a representation change.

**Do:** carry the measured threshold as a load-bearing comment with a
`BENCHMARKS.md` pointer (U7); do **not** re-derive it. `ops::Coeff`/`ops::Plan`
prepare a coefficient once when the same scalar multiplies many buffers
(`multiply`, `evaluate_many`, scale-a-domain) — amortize the preparation across
the batch, not per row.

---

## 2. Characteristic-2 `O(deg)` squaring

**What.** In characteristic 2, `(Σ aᵢ xⁱ)² = Σ aᵢ² x²ⁱ` — every cross term
`2·aᵢaⱼ` vanishes. Squaring is a coefficient **bit-spread plus a per-coefficient
square**, `O(deg)`, not the `O(deg²)` convolution a general multiply pays.

**Evidence.** `gs-engine/src/poly/arithmetic.rs:434-451` already implements
`square_into` as the spread-and-square form **[repo]**. It is the hot inner step
of `square_mod`/`pow_mod` and of the repeated-squaring in the equal-degree
root split (Cantor–Zassenhaus computes `X^(q^i)` by iterated Frobenius), so the
saving compounds along the whole root-finding path.

**Do:** keep `square`/`square_into` a first-class op distinct from
`multiply(p, p)`; route `square_mod`, `pow_mod`, and the Frobenius steps in
`roots/equal_degree.rs` through it. The oracle (§[03](03-algorithms.md)) checks
`square_into(a)` byte-identical to an independent `multiply(a, a)` — the two
share no code, so a sign-of-life bug in the spread cannot pass.

---

## 3. AFFT product above `AFFT_*_CROSSOVER`; Karatsuba fills the middle

**What.** Three product tiers by operand degree: schoolbook at the bottom,
Karatsuba in the middle, additive-FFT batched product at the top. The AFFT tier
composes `butterfly-fft` (U2); it does not contain a second FFT.

**Evidence.** The AFFT batched product already exists and is crossover-gated:
`gs-engine/src/poly/afft.rs:144-353` ships `multiply_batch_truncated[_with]`,
`PolynomialProductScratch`, and `ProductStrategy{Auto,Schoolbook,Afft}` with the
`AFFT_*_CROSSOVER` thresholds **[repo]**. Only schoolbook and AFFT exist today —
there is a **gap between them** where operands are too large for schoolbook's
`O(n²)` but too small to amortize the basis change + transform of the AFFT path.
Karatsuba–Ofman closes it at `O(n^1.585)` with no transform setup **[pub]**
(Karatsuba–Ofman 1962). This is the one *new* tier; the crossovers on both its
sides are measured, not guessed (U7).

**Do:** add `poly/karatsuba.rs` as the middle tier; wire its two crossovers into
`cost.rs` alongside the existing `AFFT_*_CROSSOVER`. Phase 5 acceptance
(§[04](04-roadmap.md)) requires all three tiers byte-identical on the same
inputs, with the crossovers recorded in `BENCHMARKS.md` and nowhere else. AFFT
is behind the default-on `fft` feature; a `--no-default-features` build keeps
schoolbook + Karatsuba only, and the crossover selector must degrade cleanly to
the two available tiers.

---

## 4. Newton-iteration truncated power-series inverse

**What.** `inverse mod x^t` by Newton doubling: `b_{k+1} = b_k (2 - a b_k) mod
x^{2^{k+1}}`, which in characteristic 2 is `b_{k+1} = b_k + b_k(1 + a b_k)`,
doubling the correct prefix each step. Cost `O(M(t))` against the `O(t²)` of a
linear schoolbook series inverse. Power-series division and reciprocal build on
it.

**Evidence.** von zur Gathen & Gerhard §9.1 is the reference for the doubling
inverse and its `O(M(t))` cost **[pub]**. This primitive exists **nowhere** in
the tree today (§[01](01-landscape.md) gap list) and is exactly what the
interleaved-RS key-equation and Padé paths downstream need. `M(t)` is whatever
tier §3 selects, so the series inverse inherits the product's crossovers for
free.

**Do:** implement in `poly/series.rs` on top of `multiply_truncated`. The oracle
(§[03](03-algorithms.md)) checks `a · a_inv ≡ 1 (mod x^t)` by an independent
truncated multiply, and cross-checks Newton doubling against a linear schoolbook
inverse — the two agree or one is wrong. Below a small `t` the schoolbook inverse
wins (no doubling overhead); carry that crossover as a comment (U7).

---

## 5. Subproduct-tree multipoint evaluation and interpolation

**What.** Evaluate a degree-`n` polynomial at `n` **arbitrary** points, or
interpolate through them, in `O(M(n) log n)` via the subproduct tree of
`∏(x − αᵢ)` and its remainder tree — against the `O(n²)` of per-point Horner or
`O(n²)` Newton. Structured (subspace/coset) domains skip the tree entirely and
compose `butterfly-fft` (U2).

**Evidence.** von zur Gathen & Gerhard §10 and Borodin–Moenck 1974 are the
reference for the subproduct-tree eval/interpolation and its `O(M(n) log n)`
cost **[pub]**. The small-`n` boundary is already measured:
`MODULE_INTERPOLATION_CROSSOVER = 8` **[repo]** — below it, Horner eval and
incremental Newton interpolation (lifted from
`gs-engine/src/interpolation/plan.rs` `build_newton_basis`,
`module.rs:470-484` `interpolate_newton_into`) win outright. The arbitrary-point
tree is the **gap**; the structured path is composition, not new code.

**Do:** build the subproduct/remainder tree in `eval/multipoint.rs`, keep
Horner/`evaluate_many` (`arithmetic.rs:138-160`) for small point sets, and route
subspace/coset domains through `TransformPlan::forward/inverse` +
`monomial↔novel` in `eval/transform.rs` under `#[cfg(feature="fft")]`. The oracle
requires the tree result equal per-point Horner, and equal `butterfly-fft`
`forward` on a subspace — three independent evaluators agreeing. This is also
where the crate **unifies the three duplicate Lagranges** it inherits (gfm's
Vandermonde synthetic-division Lagrange, `srs`'s exponent-domain Lagrange, and
gs-engine's Newton path) into one arbitrary-point primitive.

---

## 6. Root-backend selection: Chien vs equal-degree vs lifting

**What.** Three root finders, one selector. Classical **Chien search** (domain
scan, `O(n·deg)`) for the small locators of bounded-distance decoding;
`gcd(p, X^{|F|}+X)` + trace/Cantor–Zassenhaus splitting for larger-degree factor
extraction; Roth–Ruckenstein / Alekhnovich power-series lifting for the
list-decoding regime. `cost.rs` picks by degree and domain size.

**Evidence.** The equal-degree route
(`gs-engine/src/roots/field_roots.rs:104-330`) and both lifting backends
(`roots/{roth_ruckenstein,alekhnovich}.rs`, `DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER
= 20000`) already ship with selectors in `cost.rs` **[repo]**. Chien search is
the **gap** — gs-engine has always used gcd+split. Chien 1964 is the classical
`O(n·deg)` scan **[pub]**, and it wins precisely where the locator degree is a
handful and the domain is the code length: there the gcd+split's repeated
squaring `X^{q^i}` is pure overhead **[inf]**. They win in different regimes;
neither replaces the other (settled decision #4).

**Do:** add `roots/chien.rs`, wire its crossover into the existing `cost.rs`
selector, and make Phase 3 acceptance require Chien and equal-degree to agree on
root sets *and* counts, with `cost.rs` picking the cheaper on measured sizes. Do
**not** delete gs-engine's gcd+split in favour of Chien or vice versa.

---

## 7. Batch inversion for interpolation and Forney denominators

**What.** Invert a whole vector of field elements — Lagrange/Newton
denominators, Forney error-value denominators — with **one** field inversion and
`~3n` multiplies (Montgomery's product-tree trick), instead of `n` inversions.
Since `fgf::field::inv` is total (`inv(0)==0`), the batch form must preserve that
convention element-wise (U6).

**Evidence.** The pattern already appears twice in the stack —
`srs/src/tower/cauchy.rs` and `gfm/src/structured/batch_inv.rs` **[repo]** — as
the denominator-inversion step of a structured solve. It is a **field-vector**
operation, not a polynomial operation: it multiplies and inverts field elements,
touches no coefficient degree structure.

**Do:** compose whatever `fgf` exposes. If `fgf` has no batch-invert helper,
keep a thin `univariate` helper as a stopgap and **flag it upstream to `fgf`**
(it belongs in the field layer — see [08-risks.md](08-risks.md) D3), never grow
it into a second field-arithmetic home inside `univariate`. Guard the zero case
explicitly against the `inv(0)==0` oracle (U6): a batch invert that divides by a
prefix product silently breaks on a zero element.

---

## Deferred and rejected — weighed, and not to be re-proposed

| # | Proposal | Verdict | Evidence |
| - | -------- | ------- | -------- |
| R1 | Half-GCD / fast EEA (Brent–Gustavson–Yun) for gcd and truncated EEA | **Deferred, not rejected.** Its `O(M(n) log n)` beats schoolbook EEA only at degrees ≳ thousands; the constant-factor schoolbook/BM EEA wins at RS/BCH sizes (degrees in the hundreds). gfm's own plan reached the same conclusion. Gate on a consumer at degree ≳ thousands. | **[pub]/[inf]** vzGG §11; `gfm/.plans/03-algorithms.md:296-323` |
| R2 | A full GF(2^16) nibble-table bank for the eval/product lanes | **Negative.** ~9 MiB, thrashes cache; `fgf` derives its tables per call for exactly this reason. Not the crate's to build — U1 puts tables in `fgf`. | **[repo]** `fgf/src/kernel/tables.rs` |
| R3 | Frobenius-based fast root-finding beyond gs-engine's split | **Deferred.** No consumer needs it past what equal-degree already ships; adds a fourth root backend for no measured regime. | **[inf]** §[00](00-charter.md) deferred list |
| R4 | A second additive FFT inside `univariate` for structured domains | **Rejected by U2.** The transform is `butterfly-fft`'s object; the crate composes `TransformPlan` + `monomial↔novel`, never a parallel implementation. | **[repo]** ground rule 2; `butterfly-fft/src/core/transform.rs:307-680` |
| R5 | `rayon` inside the ring/root/eval hot loops | **Deferred.** Feature-gated, default off, symbol/batch axis only — never a per-coefficient or per-root task. | **[inf]** §[00](00-charter.md); §[04](04-roadmap.md) Phase 7 |

---

## Ranked: payoff ÷ cost

| # | Optimization | Payoff | Cost | Phase |
| - | ------------ | ------ | ---- | ----- |
| 1 | Packed-kernel dispatch above the measured crossover | High — U1 correctness *and* the whole ring's throughput | Very low — reuse `arithmetic.rs:573` | 1 |
| 2 | Char-2 `O(deg)` squaring | High — `square_mod`/`pow_mod`/Frobenius all ride it | Very low — exists | 1 |
| 3 | Newton-iteration series inverse `O(M(t))` vs `O(t²)` | High — the gap primitive the key-equation path needs | Medium (new) | 2 |
| 4 | Root-backend selection (Chien vs equal-degree vs lifting) | High — Chien is the cheap path for bounded-distance locators | Medium — Chien is new, selector exists | 3 |
| 5 | Subproduct-tree multipoint eval/interp `O(M(n) log n)` | High only above `MODULE_INTERPOLATION_CROSSOVER = 8` | Medium (new arbitrary-point tree) | 4 |
| 6 | AFFT product tier + Karatsuba middle tier | High — deletes `O(n²)` for large operands | Medium — AFFT exists, Karatsuba new | 5 |
| 7 | Batch inversion for interpolation / Forney denominators | Medium — `n` inversions → one | Low — compose `fgf`, else thin helper | 4 |
| R1 | Half-GCD / fast EEA | High only above degree ≳ thousands | High | deferred, gated on a consumer |

**One sentence:** the field multiplies are already `fgf`'s and already at the
ceiling, so `univariate` spends its effort on algorithms that issue *fewer*
field operations, on the char-2 shortcuts a polynomial ring uniquely has, and on
reusing the crossovers gs-engine already measured — never on re-implementing a
kernel or a transform it does not own.

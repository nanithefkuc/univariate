# Algorithms

Every algorithm the crate owns, with the **independent oracle** that decides its
correctness. An implementation is never its own oracle; where the obvious oracle
is a slower version of the same algorithm, a second, structurally different
check is named. The ring is the object (`Polynomial<F: FieldKernels>`, packed
`Vec<u8>` LE, canonical); the additive transform belongs to `butterfly-fft` and
is composed, never reimplemented (U2).

Citations are primary. Licenses are flagged where code exists: **safe** =
MIT/BSD/Apache and may be mirrored in spirit; **caution** = GPL/LGPL, read the
paper and never transcribe the code; **spec** = standards text, freely
implementable. Textbook and paper sources with no shipped code are read for the
mathematics only.

Two enduring correctness oracles cut across every section and are named once
here:

- **python `galois`** (MIT, **safe**) — `Poly` arithmetic, `.roots()`,
  `berlekamp_massey`; an API-shape and value oracle regenerated to fixed-seed
  vectors, never linked. <https://github.com/mhostetter/galois>
- **NTL** `GF2EX` / `MinPolySeq` / factoring (LGPL, **caution**) — read the
  paper, never the code; used only as a differential over FFI where a case needs
  a second independent witness.

---

## 1. Product tiers — schoolbook, Karatsuba, AFFT

One `multiply` with three interchangeable engines selected by degree
(`cost.rs`, U7). All three compute the same coefficient convolution over
GF(2^m); they differ only in cost.

**Schoolbook.** The `O(d_a·d_b)` convolution, dispatched to `fgf`'s packed
`mul_add` (AXPY) above the measured lane-bytes crossover and to a scalar
`fgf::field::Elem` loop below it — the crossover is read from gs-engine, not
re-derived (`gs-engine/src/poly/arithmetic.rs:573` `use_packed_kernel`). Extracts
`gs-engine/src/poly/arithmetic.rs:12-263`. In characteristic 2 the squaring
special case is a coefficient bit-spread, `O(deg)`, not a multiply
(`arithmetic.rs:434-451` `square_into`).

**Karatsuba.** The missing middle tier (GAP — only schoolbook and AFFT exist
today). Split each operand at the half-degree, form three products of half-size
operands, recombine with two additions; recurse to the schoolbook base case.
`O(d^{log2 3})`. In characteristic 2 the Karatsuba subtractions are XORs, so the
usual sign bookkeeping vanishes.

**AFFT.** For large operands, transform to the novel basis, pointwise-multiply
the evaluations over an additive subspace, and transform back — this *is* the
additive FFT and it is `butterfly-fft`'s (U2). Extracts the batched product
`gs-engine/src/poly/afft.rs:144-353` (`multiply_batch_truncated[_with]`,
`PolynomialProductScratch`, `ProductStrategy{Auto,Schoolbook,Afft}`), gated by
the `fft` feature and composing `TransformPlan::forward/inverse` plus
`monomial_to_novel`/`novel_to_monomial`. `O(M(n))` with `M(n) = O(n log n)`.

> Karatsuba, A., Ofman, Yu. (1962). *Multiplication of many-digit numbers by
> automatic computers.* Doklady Akad. Nauk SSSR 145, 293–294. English:
> Soviet Physics Doklady 7, 595–596 — the three-multiplication split. No shipped
> code; read the method.
>
> von zur Gathen, J., Gerhard, J. (2013). *Modern Computer Algebra*, 3rd ed.
> Cambridge University Press. DOI <https://doi.org/10.1017/CBO9781139856065> —
> §8 fast multiplication and the schoolbook↔Karatsuba↔FFT crossover analysis.
> Textbook (copyright): read, do not transcribe.
>
> The AFFT layout and kernel are `butterfly-fft`'s
> (`butterfly-fft/src/core/transform.rs:307-680`,
> `basis/convert.rs:42-183`); `univariate` composes them and owns no transform.

**Oracles.** The three tiers are each other's oracles: Karatsuba and schoolbook
must be **byte-identical** on every operand pair, and AFFT (padded to a subspace
size) must agree with both on the overlapping degrees. The base oracle beneath
all three is a naive per-coefficient convolution written once in `tests/` over
`fgf::field::Elem`, structurally unrelated to every dispatched path. Cross-check
against `galois` `Poly` products on fixed-seed vectors.

---

## 2. Euclidean division — `div_rem`

`div_rem(a, b) -> (q, r)` with `a = q·b + r` and `deg r < deg b`, monic-free
(leading coefficient inverted via `fgf::field::Elem::inv`, which is total —
`inv(0)==0`, so the zero-divisor guard is explicit, not inferred, U6). Extracts
`gs-engine/src/poly/arithmetic.rs:203-347`, which also carries `exact_divide`,
`remainder`, `monic`, `divide_by_x_power(_into)`/`x_valuation`, and the modular
wrappers `multiply_mod`/`square_mod`/`pow_mod`. Every hot form has an `_into`
scratch-owning variant (U5). Dividing by zero returns
`PolynomialError::DivisionByZero`; a non-exact `exact_divide` returns
`PolynomialError::NonExactDivision` (U9 — public boundary validates, kernels
`debug_assert!`).

> von zur Gathen & Gerhard, §2.4 (division with remainder over a field) and §9
> (the Newton-inversion route to fast division, deferred until a consumer needs
> it). DOI above.

**Oracles.** The defining identity, checked by the *independent* multiply/add of
§1: `q·b + r == a` and `deg r < deg b`, on randomized `(a, b)` at every degree
including `deg a < deg b` (quotient zero) and `b` a unit. The multiply used for
the check is the naive convolution oracle, never `div_rem`'s own back-multiply.

---

## 3. GCD and extended GCD with Bézout cofactors

Two routines share one loop. The **plain gcd** already exists — monic Euclidean
remainder sequence, `gs-engine/src/poly/arithmetic.rs:271-281` — and is
extracted unchanged. The **extended gcd** is a GAP (only the plain gcd ships
today): it accumulates the Bézout cofactors `s, t` alongside the remainder
sequence so that `s·a + t·b == g`, returning the monic `g` and both cofactors.
The accumulation is the same 2×2 unimodular transform the truncated EEA of §4
stops early; extended gcd runs it to `g` (remainder zero at the next step).

> von zur Gathen & Gerhard, §3.2 (the extended Euclidean algorithm) and §6
> (Bézout relation, cofactor degree bounds `deg s < deg b − deg g`,
> `deg t < deg a − deg g`). DOI above. Read for the degree bounds that make the
> cofactors unique.

**Oracles.** The Bézout identity `s·a + t·b == g` verified by the §1 multiply
and §2 add; divisibility `g | a` and `g | b` verified by `div_rem` returning
zero remainder — a different code path from the one that produced `g`. The
cofactor degree bounds are asserted exactly (not as predicates). `gcd(a, 0) == a`
made monic, and `gcd(a, a) == a` monic, are fixed cases. Cross-check `g` against
`galois` `gcd` on fixed vectors.

---

## 4. Truncated / partial EEA — the key-equation primitive

The reusable Padé / connection-polynomial primitive (GAP — exists nowhere in the
tree). Run the extended Euclidean algorithm on `(x^t, S(x))` — or on `(a, b)` —
accumulating the 2×2 unimodular transform, and **stop the instant the remainder
degree drops below `t`**. The stopped remainder and its cofactor are the Padé
approximant / the error-locator–evaluator pair of the key equation. Dornstetter
proved this stopped-EEA is exactly Berlekamp–Massey on the same sequence, so the
one routine serves both the Euclidean and the shift-register views; the actual
BM/decoding *orchestration* is `syndrome-engine`'s, and this crate supplies only
the polynomial primitive (settled decision #2; the gfm Hankel/BM seam is an open
cross-plan reconciliation, [`08-risks.md`](08-risks.md) D1 — `univariate` owns
the F[x] EEA, gfm re-scopes its unshipped Hankel claim to the scalar-sequence /
Wiedemann view or composes this).

> Sugiyama, Y., Kasahara, M., Hirasawa, S., Namekawa, T. (1975). *A method for
> solving the key equation for decoding Goppa codes.* Information and Control
> 27(1), 87–99. DOI <https://doi.org/10.1016/S0019-9958(75)90090-X> — the
> Euclidean key-equation solver.
>
> Massey, J. L. (1969). *Shift-register synthesis and BCH decoding.* IEEE Trans.
> Inf. Theory 15(1), 122–127. DOI <https://doi.org/10.1109/TIT.1969.1054260> —
> the LFSR-synthesis view the stopped EEA is equivalent to.
>
> Dornstetter, J.-L. (1987). *On the equivalence between Berlekamp's and
> Euclid's algorithms.* IEEE Trans. Inf. Theory 33(3), 428–431. DOI
> <https://doi.org/10.1109/TIT.1987.1057299> — the equivalence that lets one
> routine serve both.
>
> Brent, R. P., Gustavson, F. G., Yun, D. Y. Y. (1980). *Fast solution of
> Toeplitz systems of equations and computation of Padé approximants.* J.
> Algorithms 1(3), 259–295. DOI <https://doi.org/10.1016/0196-6774(80)90013-9> —
> the O(n log² n) half-GCD version. **Deliberately deferred**, see §10; recorded
> so the reason survives.

**Oracles.** Two, structurally independent:
1. **Padé identity.** The reconstructed rational function `t(x)/s(x)` matches the
   input series `S(x)` modulo `x^t` — checked by the truncated multiply of §9,
   not by rerunning the EEA.
2. **Dornstetter equivalence.** A reference Berlekamp–Massey, written once in
   `tests/` over `fgf::field::Elem` and never optimized, run on the same
   sequence must return the same connection polynomial (up to the documented
   normalization). Two algorithms, one answer.

Plus the LFSR property directly: the returned connection polynomial reproduces
the sequence. Cross-check against `galois::berlekamp_massey` on fixed vectors.

---

## 5. Evaluation — Horner and subproduct-tree multipoint

**Horner.** `evaluate(p, α)` and `evaluate_many(p, &[α])` are the single-point
and small-set path, extracted from `gs-engine/src/poly/arithmetic.rs:138-160`.
Scalar `fgf::field::Elem` fused multiply-add per coefficient; the batch form
lifts to `fgf`'s packed `mul_elementwise` when the point set is wide enough
(U1). The Hasse-derivative evaluation `evaluate_hasse` and `hasse_derivative`
travel with it (`arithmetic.rs:434-571`, feeding `hasse`).

**Subproduct-tree multipoint.** For a large *arbitrary* point set, build the
subproduct tree of `∏(x − αᵢ)`, then recursively reduce `p` modulo each subtree
node; the leaves are the evaluations. `O(M(n) log n)` versus `O(n·deg)` Horner,
with the crossover measured (small-`n` stays Horner). This arbitrary-point path
is the GAP; the structured-domain path composes `butterfly-fft`
(`TransformPlan::forward` = multipoint evaluation over an additive subspace, U2)
and owns no transform.

> von zur Gathen & Gerhard, §10.1 (subproduct tree, fast multipoint evaluation).
> DOI above.
>
> Borodin, A., Moenck, R. (1974). *Fast modular transforms.* J. Computer and
> System Sciences (evaluation/interpolation by the subproduct tree). DOI
> <https://doi.org/10.1145/321662.321664> — the tree construction.
>
> Structured-subspace evaluation is `butterfly-fft::TransformPlan::forward`
> (`butterfly-fft/src/core/transform.rs:307-680`).

**Oracles.** The subproduct-tree result equals per-point Horner on every point,
Horner being the structurally simpler and separately tested routine. Over an
additive subspace, the multipoint result equals `butterfly-fft`'s `forward`
after `monomial_to_novel` — a completely different (transform) code path. Both
oracles run on fixed-seed point sets.

---

## 6. Interpolation — Newton, Lagrange, inverse-FFT

Recover the minimal-degree `p` (`deg p < n`) from `n` point/value pairs, by three
routes that must agree.

**Newton.** Incremental Newton-basis interpolation with divided differences,
lifted from `gs-engine/src/interpolation/plan.rs` (`build_newton_basis`, `O(n²)`)
and `module.rs:470-484` (`interpolate_newton_into`). The *univariate* parts move
here; the poly-**matrix** / weak-Popov module interpolation stays in gs-engine
and keeps calling `gfm::weak_popov` (item 6, LIFT). Small-`n` default
(`MODULE_INTERPOLATION_CROSSOVER=8`).

**Lagrange.** Arbitrary points via the subproduct tree of §5: build `M = ∏(x −
αᵢ)`, evaluate `M'` at each `αᵢ`, and combine. `O(M(n) log n)`. This is the one
canonical Lagrange; it retires the three duplicates the ecosystem grew — gfm's
Vandermonde hand-rolled synthetic-division Lagrange
(`gfm/src/structured/vandermonde.rs:127-160`, migrate) and srs's exponent-domain
Λ(x)/((x⊕d)·Λ'(d)) (`srs/src/afft/generator.rs:11-45`, noted not extracted).

**Inverse-FFT.** Over an additive subspace, interpolation *is*
`butterfly-fft::TransformPlan::inverse` followed by `novel_to_monomial` (U2). The
`fft` feature gates it.

> von zur Gathen & Gerhard, §5.3 (Newton form, divided differences) and §10.2
> (fast interpolation by the subproduct tree). DOI above.
>
> Borodin–Moenck (1974), DOI above — the interpolation companion to §5.
>
> Inverse subspace interpolation is `butterfly-fft::TransformPlan::inverse`
> (`butterfly-fft/src/core/transform.rs:307-680`).

**Oracles.** Round-trip: interpolate, then `evaluate_many` (§5) at the original
points returns the input values exactly, and `deg p < n`. The three routes must
agree with each other on shared domains — Newton vs Lagrange on arbitrary
points, all three on a subspace — which is the differential the ecosystem's
three duplicate Lagranges never had. Cross-check against `galois` interpolation
on fixed vectors.

---

## 7. Chien search

The classical root search by domain scan (GAP — gs-engine uses gcd+split
instead). Evaluate the locator at every element of the field or a coset by an
incremental update: stepping `α → α·γ` multiplies the `j`-th coefficient's
running term by `γ^j`, so each successive evaluation costs one packed
`mul_elementwise` plus a reduction rather than a fresh Horner. `O(|domain|·deg)`.
This is the cheap path for the *small* locators of bounded-distance decoding;
`cost.rs` selects it against equal-degree factorization (§8), which wins for
larger-degree factor extraction — both ship, neither replaces the other (settled
decision #4). Enumeration order is a frozen wire property (U8): consumers map
roots to positions.

> Chien, R. T. (1964). *Cyclic decoding procedures for Bose–Chaudhuri–
> Hocquenghem codes.* IEEE Trans. Inf. Theory 10(4), 357–363. DOI
> <https://doi.org/10.1109/TIT.1964.1053699> — the root-search-by-scan and its
> incremental update. **spec**-adjacent classical method; freely implementable.

**Oracles.** Every enumerated root `r` satisfies `evaluate(p, r) == 0` by the
independent Horner of §5. The root *set* and its *count* match §8's equal-degree
factorization on the same polynomial — a structurally unrelated algorithm — and
the count respects the degree bound. A companion **linearized / affine
(2-linearized) root solver** (GAP) handles the affine-polynomial special case in
closed form; its roots are verified by the same `evaluate == 0` oracle and must
match the Chien enumeration on the affine inputs.

---

## 8. Root-finding — equal-degree factorization and power-series lifting

**Equal-degree (base-field roots).** `base_field_roots` computes `gcd(p, X^{|F|}
+ X)` to isolate the degree-1 factors, then Cantor–Zassenhaus trace/Frobenius
splitting recovers the individual roots, reusing the scratch-owning `mul_mod`/
`square_mod`. Extracted verbatim from `gs-engine/src/roots/field_roots.rs:104-330`.

**Roth–Ruckenstein / Alekhnovich lifting.** Power-series root lifting for the
list-decoding regime — lift a root of a bivariate/`x`-adic reduction coefficient
by coefficient. Both backends extract from
`gs-engine/src/roots/{roth_ruckenstein,alekhnovich}.rs` with their
`AlekhnovichScratch` and the measured crossover
`DEFAULT_ROTH_RUCKENSTEIN_CROSSOVER=20000` (`cost.rs`, U7).

> Cantor, D. G., Zassenhaus, H. (1981). *A new algorithm for factoring
> polynomials over finite fields.* Math. Comp. 36(154), 587–592. DOI
> <https://doi.org/10.2307/2007663> — the gcd(p, X^q − X) + splitting route.
>
> Berlekamp, E. R. (1968, rev. 2015). *Algebraic Coding Theory.* World
> Scientific — equal-degree factorization and root-finding over GF(2^m).
>
> Roth, R. M., Ruckenstein, G. (2000). *Efficient decoding of Reed–Solomon codes
> beyond half the minimum distance.* IEEE Trans. Inf. Theory 46(1), 246–257. DOI
> <https://doi.org/10.1109/18.817522>
>
> Alekhnovich, M. (2005). *Linear Diophantine equations over polynomials and
> soft decoding of Reed–Solomon codes.* IEEE Trans. Inf. Theory 51(7),
> 2257–2265. DOI <https://doi.org/10.1109/TIT.2005.850102>

**Oracles.** Every returned root verified by the independent `evaluate` of §5
(`evaluate(p, r) == 0`). The two lift backends cross-check against each other —
Roth–Ruckenstein and Alekhnovich must return the same root set on the same input
(already gs-engine practice). Equal-degree root sets and counts agree with the
Chien scan of §7 on the base-field case. Cross-check against `galois` `.roots()`
on fixed vectors; NTL factoring over FFI as a second witness where a case needs
one.

---

## 9. Truncated power-series inverse

`inverse_mod_x_t(a, t)` returns `a⁻¹ mod x^t` for `a` with a nonzero constant
term (GAP — no truncated-series algebra exists in the tree today). Newton
doubling: seed with the inverse of the constant term (`fgf::field::Elem::inv`,
`inv(0)==0` guarded explicitly, U6), then iterate `b ← b·(2 − a·b) mod x^{2k}`,
doubling the correct prefix each step. In characteristic 2 the `2` vanishes, so
the step is `b ← b·(1 − a·b) = b + a·b²` truncated — an even cheaper update.
`O(M(t))` versus `O(t²)` for the schoolbook long-division inverse. Series
division and the reciprocal build on it; `mul_trunc` already exists as
`multiply_truncated` (§1, `arithmetic.rs`). This is the primitive
`syndrome-engine` and the interleaved-RS path in `contort` need.

> von zur Gathen & Gerhard, §9.1 (Newton iteration for power-series inversion,
> the doubling schedule and its cost). DOI above.

**Oracles.** The defining identity `a · a_inv ≡ 1 (mod x^t)` checked by the
truncated multiply of §1 — not by the inverse's own internal product. The Newton
doubling result must equal a **linear schoolbook series inverse** written once in
`tests/` (solve for one coefficient at a time), a structurally different
algorithm giving the same series. `inverse_mod_x_t` of a unit is that unit's
reciprocal; a zero constant term returns the documented error, not a silent zero.

---

## 10. What we deliberately do not implement

Recorded with the reason, so it is not relitigated.

**Half-GCD / fast EEA (Brent–Gustavson–Yun).** The O(n log² n) recursive
truncated-EEA. At RS/BCH sizes — degrees in the hundreds — the constant-factor
schoolbook / BM EEA of §4 wins outright, and gfm's own plan reaches the identical
conclusion (`gfm/.plans/03-algorithms.md:296-323`). Gated on a consumer at degree
≳ thousands, not before.
> Brent, Gustavson, Yun (1980), DOI <https://doi.org/10.1016/0196-6774(80)90013-9>.

**Frobenius-based fast root-finding beyond gs-engine's route.** The
equal-degree + lifting backends of §8 already cover the list-decoding regime;
faster Frobenius/Kaltofen–Shoup factorization buys nothing at these degrees and
adds a large implementation surface. Gated on a consumer that profiles into it.

**The additive FFT itself, novel/Cantor bases, transform buffers.** Not ours at
all — `butterfly-fft` owns the layout and kernel (U2). A plan that claimed to
accelerate subspace evaluation by writing a second additive FFT would be wrong.

**Bivariate polynomials and the GS interpolation/list-decoding machinery.** Stay
in gs-engine (`gs-engine/src/poly/bivariate.rs`); `univariate` is strictly
univariate and the bivariate layer is a consumer of it (settled decision #5).

**Matrices of polynomials / weak-Popov row reduction.** gfm's object
(`gfm/src/poly.rs`); `univariate` supplies the F[x] *element* ops gfm composes,
not the reduction.

---

## Reading order for the implementer

1. **von zur Gathen & Gerhard, §2 and §9** — division with remainder, then the
   Newton-inversion machinery that §9 (series inverse) and later fast division
   share. Read before writing `divide.rs` and `series.rs`.
2. **Karatsuba–Ofman 1962** — the three-multiplication split, then the §1
   crossover analysis in vzGG §8. Build the middle tier against schoolbook.
3. **von zur Gathen & Gerhard, §10** + **Borodin–Moenck 1974** — the subproduct
   tree unifies multipoint evaluation (§5) and interpolation (§6). One data
   structure, two directions.
4. **Sugiyama 1975** + **Dornstetter 1987** + **Massey 1969** — build the
   truncated EEA (§4) once; the Dornstetter equivalence is the oracle and the
   reason `syndrome-engine` need not grow its own.
5. **Chien 1964** — the incremental scan (§7), and where `cost.rs` picks it over
   equal-degree.
6. **Cantor–Zassenhaus 1981** + **Berlekamp 1968** — the equal-degree route
   (§8); then **Roth–Ruckenstein 2000** and **Alekhnovich 2005** for the lifting
   backends, read as the pair that cross-checks the other.
7. **`butterfly-fft/src/core/transform.rs` and `basis/convert.rs`** — the
   `TransformPlan` + `monomial↔novel` surface `univariate` composes for §1 (AFFT
   tier), §5 (subspace eval), and §6 (inverse-FFT interpolation).

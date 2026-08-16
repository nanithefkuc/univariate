# Baselines and reference sources

Two questions. **What do we measure against**, and **what do we implement from**.
They have different answers, and the difference is mostly licensing. The
algorithm bibliography lives in [`03-algorithms.md`](03-algorithms.md) and
[common-dossier]; this document indexes it and states the license posture.

## License posture

This crate is **MIT**. Algorithms are not copyrightable; implementations are.

| Class | Projects | Rule |
| --- | --- | --- |
| **Safe** — may mirror closely | python `galois` (MIT), `reed-solomon-erasure` (MIT), `reed-solomon-simd` (MIT AND BSD-3), `klauspost/reedsolomon` (MIT) | Read the code; mirror the structure; attribute. |
| **Caution** — read the paper, never the code | NTL `GF2EX` (LGPL-2.1+), FLINT `fq_nmod` (LGPL-3.0+), SageMath (GPLv2+), PARI/GP (GPLv2+) | Benchmark against them. Cite the papers. Do not transcribe code or comments. |
| **Spec** | IEEE/textbook decode pipelines (Chien, Forney, Euclidean key equation) | Implement from the text. |
| **Cite only** | Magma (proprietary) | Published numbers as a reference point. |
| **In-repo** | `gs-engine/src/{poly,roots,interpolation}`, `gs-engine/src/{domain,cost}.rs` | The extraction source — this *is* the reference implementation; the differential is its own existing test suite. |

The point of the table: the mature polynomial ring `univariate` extracts is
**already MIT and already ours** (`gs-engine`), so the primary correctness
oracle is a differential against the code being lifted, not a foreign library.
The external baselines exist to check the **gap primitives** — cofactor EEA,
truncated EEA, Chien, power-series inverse — that have no in-repo predecessor.

---

## Part 1 — What we benchmark against

### The three primary baselines

**1. python `galois` (MIT) — the ergonomics and correctness comparator.**

<https://github.com/mhostetter/galois> — MIT, actively maintained, same field
range (all GF(p^m), so every `Gf8B/Gf8D/Gf16/Gf32/Gf64` and Fan–Paar level has a
counterpart), and a NumPy-shaped API that is the closest thing to a *reference
spec* for the operations `univariate` owns:

- `galois.Poly` arithmetic (`+`, `*`, `//`, `%`, `divmod`, `**`), `.derivative()`
- `Poly.roots()` — the root-finding oracle for Chien / equal-degree / lifting
- `galois.berlekamp_massey` — the independent oracle for the truncated-EEA key
  solver via the Dornstetter equivalence (§[03](03-algorithms.md))
- `Poly.reverse()`, modular arithmetic, and Lagrange interpolation helpers

It is textbook and Numba-JIT'd — expect it slower than the extracted Rust — so it
is an **API-shape and correctness** reference, not a performance target. Its MIT
license makes it the one baseline we may read closely and mirror in spirit.

**2. NTL `GF2EX` (LGPL, caution) — the independent large-field oracle and floor.**

<https://libntl.org> — `GF2EX` is dense univariate polynomials over GF(2^m) with
the full ring: `mul`, `DivRem`, `GCD`, `XGCD` (extended Euclid with cofactors —
the direct oracle for the crate's *new* cofactor EEA), `MinPolySeq` (the BM /
key-equation oracle), `FindRoots`/factoring, and `eval`. It reaches the wide
fields galois exercises slowly. **Caution:** LGPL — benchmark against it and cite
the papers; do **not** transcribe its code or comments. "If the crate's XGCD or
root set disagrees with NTL's on the same input, the crate is wrong."

**3. The crate's own naive tier — the differential that writes itself.**

Every fast path ships beside a naive reference in the same crate: schoolbook
convolution behind Karatsuba/AFFT, per-point Horner behind the subproduct tree,
linear schoolbook series-inverse behind Newton doubling, gcd+split behind Chien.
The naive tier is not dead weight — it **is** the primary oracle (§[03](03-algorithms.md)):
each fast algorithm is required byte-identical to the naive one it replaces, and
the two share no code. This is cheaper and stronger than any foreign library
because it runs over the exact same `Polynomial<F>` representation and the exact
same field, so a disagreement localizes to one algorithm, not to a
representation or field-model mismatch.

### Secondary and cross-check baselines

| Project | Use | Note |
| --- | --- | --- |
| **`gs-engine`** (in-repo) | The extraction source and the reference consumer. Its existing `poly`/`roots`/`interpolation` tests are the differential for Phase 1–6; Phase 6's cutover requires its suite to pass unchanged against `univariate`. | MIT, ours. The strongest oracle we have — same code, pre- and post-extraction. |
| **FLINT** `fq_nmod_poly` <https://flintlib.org> | Cross-check for the wide fields and for fast division / half-GCD asymptotics, should a consumer ever push degrees into the thousands (deferred, [06](06-optimizations.md) R1). | LGPL-3.0+, caution — benchmark and cite only. |
| **SageMath** | Zero-build-effort driver over GF(2^m)`[x]` (`R.<x> = GF(2^8)[]`) for spot correctness of gcd, xgcd, roots, factor, interpolate. Python overhead is nonzero; use `cputime()`. | GPLv2+, caution. |
| **PARI/GP** | Scriptable checker (`gp -q`): `gcdext`, `polroots`/`factorff`, `polinterpolate` over `ffgen`. Weak as a speed rival, convenient as a second oracle. | GPLv2+, caution. |
| **`reed-solomon-erasure`** <https://github.com/rust-rse/reed-solomon-erasure> | Perf floor and correctness datapoint for the GF(2^8) evaluation/Horner path — the only Rust crate shipping a real GF matrix inverse, useful as a "what users have today" marker. | MIT. **Unmaintained** — "Looking for new owners". |
| **`reed-solomon-simd`** <https://github.com/AndersTrier/reed-solomon-simd> | Bulk-throughput sanity bar for the transform product tier: Leopard-style `O(n log n)` GF(2^16) FFT, so it bounds what the `butterfly-fft`-composed AFFT product should cost. | MIT AND BSD-3. Erasure-only — not a decode-pipeline comparison. |
| **`klauspost/reedsolomon`** <https://github.com/klauspost/reedsolomon> | Bulk GF(2^8)/GF(2^16) row-operation throughput reference, to catch an evaluation/product kernel leaving SIMD on the floor. | MIT. Erasure-only. |
| **Magma** | Cite only. Published `GF2EX` gcd/factor/interpolate numbers as a reference point for the wide fields. | Proprietary. |

### Is this greenfield?

**Partly — and the split is the whole point.**

- **The ring is not greenfield; it is trapped.** A full, `fgf`-dispatched,
  packed-byte univariate library already exists in this repo — `gs-engine/src/poly/`
  (`Polynomial<F>` + `arithmetic.rs`) and `gs-engine/src/roots/`. It cannot be
  reused by `syndrome-engine`, `contort`, `funcfield`, `hasse`, or `reed-muller`
  because it is a private module of the Guruswami–Sudan engine. `univariate` is
  the extraction that turns it into an L1 object. This is the mirror image of
  `gfm`, which collapsed six copies; here there is one canonical copy and several
  thin duplicates (gfm's Vandermonde synthetic-division Lagrange, `srs`'s
  exponent-domain Lagrange) to retire.
- **The classical decoder primitives *are* greenfield, in Rust and mostly in the
  permissive world at large.** Extended Euclid with Bézout cofactors, the
  truncated/partial EEA that solves a key equation, classical Chien search, and
  truncated power-series inversion (`inv mod x^t` / Newton / series division)
  exist **nowhere** in this tree, and the syndrome→BM/EEA→Chien→Forney path is
  rare in permissive Rust — most RS crates are erasure-only (galois is the MIT
  exception, and it is Python). These are the fresh build.

The defensible claim for `univariate` is therefore precise: a safe-Rust,
permissively licensed univariate polynomial ring over GF(2^m) that composes `fgf`
for arithmetic and `butterfly-fft` for structured transforms, lifting the
one mature in-repo implementation to L1 and adding the classical
key-equation/root primitives that no permissive Rust crate provides.

---

## Part 2 — What we implement from

Grouped by the component they inform. Full annotations, oracles, and reading
order are in [`03-algorithms.md`](03-algorithms.md); this is the index. Every
entry cites only the sanctioned bibliography (see [common-dossier]).

| Component | Primary sources | Code to read |
| --- | --- | --- |
| Ring, division, plain gcd | von zur Gathen & Gerhard §2–3 [DOI](https://doi.org/10.1017/CBO9781139856065) | `gs-engine/src/poly/arithmetic.rs` (ours, the extraction source); NTL `GF2EX` (caution) |
| Schoolbook / Karatsuba / AFFT product | Karatsuba–Ofman 1962, *Doklady* 145, 293–294; von zur Gathen & Gerhard §8 | `gs-engine/src/poly/afft.rs` (ours, AFFT tier); `butterfly-fft/src/core/transform.rs` (ours, the transform) |
| Extended Euclid + Bézout cofactors | von zur Gathen & Gerhard §3, §6 [DOI](https://doi.org/10.1017/CBO9781139856065) | galois `berlekamp_massey` / NTL `XGCD` (oracles) — **new, no in-repo predecessor** |
| Truncated / partial EEA (key equation) | Sugiyama–Kasahara–Hirasawa–Namekawa 1975 [DOI](https://doi.org/10.1016/S0019-9958(75)90090-X); Massey 1969 [DOI](https://doi.org/10.1109/TIT.1969.1054260); Dornstetter 1987 [DOI](https://doi.org/10.1109/TIT.1987.1057299) | galois `berlekamp_massey`; NTL `MinPolySeq` (oracles, caution) — **new** |
| Truncated power series (`inv mod x^t`, division) | von zur Gathen & Gerhard §9.1 (Newton inversion) [DOI](https://doi.org/10.1017/CBO9781139856065) | — **new, no in-repo predecessor** |
| Multipoint eval + interpolation (subproduct tree) | von zur Gathen & Gerhard §10; Borodin–Moenck 1974 [DOI](https://doi.org/10.1145/321662.321664) | `gs-engine/src/interpolation/plan.rs`, `module.rs:470-484` (Newton path, ours); `butterfly-fft` (structured domains, ours) |
| Chien search | Chien 1964, IEEE Trans. IT 10(4), 357–363 [DOI](https://doi.org/10.1109/TIT.1964.1053699) | galois `Poly.roots()` (oracle) — **new, gs-engine uses gcd+split** |
| Equal-degree factorization / base-field roots | Cantor–Zassenhaus 1981 [DOI](https://doi.org/10.2307/2007663); Berlekamp, *Algebraic Coding Theory* 1968 | `gs-engine/src/roots/field_roots.rs:104-330` (ours) |
| Power-series root lifting | Roth–Ruckenstein 2000 [DOI](https://doi.org/10.1109/18.817522); Alekhnovich 2005 [DOI](https://doi.org/10.1109/TIT.2005.850102) | `gs-engine/src/roots/{roth_ruckenstein,alekhnovich}.rs` (ours) |
| Evaluation domains + backend selection | von zur Gathen & Gerhard §8, §10 (crossover analysis) | `gs-engine/src/{domain,cost}.rs` (ours) |
| Forney / errors-and-erasures context (consumer) | Forney 1965 [DOI](https://doi.org/10.1109/TIT.1965.1053825); Blahut 1983/2003; Lin & Costello 2004 | consumed by `syndrome-engine`; `univariate` supplies `formal_derivative` + eval + batch-inv |

---

## Benchmark harness plan

`external-bench/` at the crate root, out of tree, following the stack's
precedent of a comparison directory that is **not** part of the crate and **not**
in CI's required path.

1. **`galois/`** — a Python script driving `galois.Poly` over each field level:
   gcd, xgcd, `roots()`, `berlekamp_massey`, interpolate, evaluate. Emits the API
   comparison and honest wall-clock, so the ergonomics claim carries a number.
2. **`ntl/`** — a small C++ shim exposing `GF2EX` `XGCD`, `MinPolySeq`,
   `DivRem`, `GCD`, `FindRoots` behind one `extern "C"` surface, gated on NTL
   being installed and **skipping loudly** when it is not. The wide-field oracle.
3. **`differential/`** — the in-crate naive tier vs the fast tier over the same
   `Polynomial<F>` and the same fixed-seed inputs. This is the primary oracle and
   the one that runs in CI; the two above are opt-in.
4. **`gs-engine/`** — Phase 6's headline: re-run gs-engine's own `poly`/`roots`
   suite against the extracted `univariate` and assert its benches land within
   5% (§[04](04-roadmap.md) Phase 6).

Every result lands in `BENCHMARKS.md` with hardware, rustc version, library
versions, and the exact command line — and nowhere else (U7).

[common-dossier]: ../../.plans/common-dossier.md

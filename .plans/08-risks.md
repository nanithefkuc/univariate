# Risks and open decisions

Risks first, then the decisions that need a human call. Each risk names its
trigger, its blast radius, and the mitigation already built into the plan. Each
decision has a default so work is never blocked on the answer.

## Risks

### R1 — The extraction diverges from `gs-engine`

**Trigger.** The mature `Polynomial<F>` + `roots/` lives today as `gs-engine`'s
private `poly/` module. Lifting it into L1 while `gs-engine` keeps shipping
creates a dual-maintenance window in which the two copies drift on
normalization, the zero-polynomial convention, the packed-vs-scalar crossover,
or the `PolynomialError` variants.

**Blast radius.** Correctness and interop: a consumer migrating from one to the
other that hits a convention difference sees a behavior change with a green test
suite. The whole point of the extraction is to have one ring; a drift undoes it.

**Mitigation, in three layers.**

1. *Single source of truth during the window.* The L1 copy is authoritative; any
   bug fix lands there first and is pulled into `gs-engine` by rev pin, never the
   reverse. The `gs-engine` private copy is treated as a vendored snapshot, not a
   fork.
2. *Convention parity at extraction.* `00-charter.md` U3 and the naming table in
   `05-conventions.md` pin the conventions to `gs-engine`'s current shape (packed
   `Vec<u8>` LE, empty buffer = zero, `normalize()` trims high zeros, `inv(0)==0`
   inherited). The lift does not retcon them.
3. *Phase 6 is a clean cutover.* `gs-engine` deletes `src/poly/` and `src/roots/`
   and imports `univariate`; its existing test suite is the differential, and the
   acceptance criterion is that it passes unchanged. No `#[deprecated]`, no
   re-export shims. After Phase 6 the dual-maintenance window is closed.

### R2 — The `gfm` Hankel / truncated-EEA seam (highest)

**Trigger.** `gfm`'s planning set claims the key-equation "one engine" as its
in-scope object — "Hankel/Toeplitz solve — `minimal_lfsr` / Berlekamp–Massey /
truncated EEA … Massey's algorithm *is* a Hankel solver; Dornstetter proved it is
the extended Euclidean algorithm" (`gfm/.plans/00-charter.md:49`), with a planned
`poly/hankel.rs` (`02-architecture.md:50`) and a full algorithm spec
(`03-algorithms.md:296-323`). This crate's settled decision #2 puts the
polynomial truncated-EEA / connection-polynomial primitive in `univariate`, as
gcd-with-cofactors over F[x] adjacent to the `div_rem`/`gcd` it already owns.

**Blast radius.** Ground rule 1: if both crates implement the primitive, the
ecosystem grows a duplicate — the exact defect `gfm` was created to remove. The
two framings are not even clearly the same object: `gfm`'s is the
linear-algebra / Wiedemann-Krylov view over a scalar sequence (oracle: a dense
`Ple` Hankel solve), while `univariate`'s is polynomial algebra (oracle: a Padé
reconstruction matching the input series mod `x^t`).

**Mitigation, in three layers.**

1. *No code collision today.* `gfm`'s `hankel.rs` is **unshipped**: the working
   tree has only `gfm/src/poly.rs` = shifted weak-Popov row reduction. The seam
   is between two *plans*, not two implementations. Reconciling before either
   writes the code is cheap.
2. *The object boundary decides it.* Truncated EEA over F[x] is a polynomial-ring
   operation; `univariate` owns the ring. `gfm`'s legitimate claim is the
   *matrix* / black-box-Wiedemann use (minimal polynomial of a sparse operator),
   which is linear algebra and stays in `gfm`. The two framings then compose
   rather than collide: `gfm`'s `minimal_lfsr` on a scalar sequence can call
   `univariate`'s truncated EEA if it wants the Dornstetter reduction, or solve
   the Hankel system directly if it wants the matrix view.
3. *Honest accounting.* This risk and D1 record the seam out loud; the crate
   does not silently override a sibling plan. The reconciliation is a graph-and-
   charter amendment to `gfm/.plans`, not a fait accompli in `univariate/src`.

### R3 — The `fft` feature leaks into the core ring

**Trigger.** `butterfly-fft` is an optional, default-on dependency behind the
`fft` feature. The core ring (add, multiply, divide, gcd, EEA, Horner eval,
Chien, power-series) must build and test `--no-default-features` with no
`butterfly-fft` in the tree, because `syndrome-engine` depends on exactly that
configuration.

**Blast radius.** A `#[cfg(feature = "fft")]` path that accidentally gate-keeps
a core API forces every transform-free consumer to take the whole FFT stack, or
to see compile errors. Quietly, with a green suite under `--all-features`.

**Mitigation.** Phase 0's `no-fft` CI job builds and tests
`--no-default-features`; Phase 0's `deps` job asserts the dependency tree at
`--no-default-features` contains exactly `fgf` and no `butterfly-fft`. The
structured-domain eval/interp is a separate module (`eval/transform.rs`) with
its own types; it is not the implementation of a core-ring trait method.

### R4 — Determinism drift in root/interpolation order

**Trigger.** Consumers map roots to positions (`syndrome-engine`'s Chien→column
map, `contort`'s interleaved-RS locator roots→column map). If `base_field_roots`
or Chien returns roots in an order that depends on the trace-splitting seed, the
Karatsuba-vs-schoolbook tie-break, or a `HashMap` iteration, the mapping breaks
silently.

**Blast radius.** Decoder interop: two decoders with the same locator produce
corrections at different positions. This is the "Determinism is a wire property"
invariant from the root `AGENTS.md`, and it is the one a refactor is most likely
to break without any test failing.

**Mitigation.** U8 makes root enumeration order a documented, frozen order
(sorted by position, or by the canonical field ordering `fgf` exposes), and
Phase 3's acceptance checks it as a property on randomized inputs across
backends. `cost.rs`'s backend selectors may choose a different *algorithm* but
must produce the same *ordered* root set. No `HashMap`/`HashSet` in a root or
interpolation return path.

---

## Open decisions

These need a human call; each has a recommendation, and each is defaulted so
work is never blocked on the answer.

### D1 — Reconcile the truncated-EEA seam with `gfm`

Two framings, one primitive. **(A) `univariate` owns the polynomial truncated
EEA**; `gfm`'s `hankel.rs` is either dropped or re-scoped to the genuine
black-box-Wiedemann minimal-polynomial routine on a scalar sequence, composing
`univariate` when it wants the Dornstetter reduction. **(B) `gfm` owns it** as a
Hankel solver and `univariate` exposes only plain + extended gcd, leaving the
key-equation primitive to the matrix view.

**Default: (A).** The primitive is gcd-with-cofactors over F[x] — a polynomial-ring
operation — and `univariate` already owns `div_rem` and `gcd` (`gs-engine/src/poly/arithmetic.rs:206,260`).
Putting it here gives every BCH/RS/interleaved-RS caller one L1 home without
pulling `gfm`'s matrix stack. `gfm`'s matrix view composes it rather than
re-implementing it, preserving one-elimination-per-domain on both sides.

*Revisit when* `gfm` demonstrates a black-box sparse-solve consumer that needs its
own minimal-polynomial routine on a scalar sequence *and* that routine is
genuinely cheaper as a direct Hankel solve than as a call into `univariate`'s
EEA. Then the reconciliation is an amendment to `gfm/.plans` first, in a
separate discussion, not a silent second implementation here.

### D2 — Is `fft` default-on?

The core ring needs no transform; `syndrome-engine` builds
`--no-default-features`. But `gs-engine`, `hasse`, `funcfield`, and
`reed-muller` all want the structured-domain path, and an L1 crate that makes
its most useful feature opt-in surprises every consumer.

**Default: yes, `fft` is in `default`.** `syndrome-engine` opts out; everyone
else gets the FFT paths by default. The `no-fft` CI job proves the opt-out
builds and the core API does not leak.

*Revisit only if* a second transform-free consumer appears wanting `fft` off by
default — at which point the question is whether *any* L1 crate should default
its optional transform on, and the answer is probably still yes for the same
reason.

### D3 — Does the field batch-invert helper move to `fgf`?

Montgomery batch inversion (one inversion plus `3(n−1)` multiplications) is a
field-vector operation, not a polynomial operation. It appears twice in the tree
today — `srs/src/tower/cauchy.rs:120,131` and `gfm/src/structured/batch_inv.rs`
— and `univariate`'s interpolation denominators and Forney denominators (via
`syndrome-engine`) need it.

**Default: compose whatever `fgf` exposes; if `fgf` has none, keep a thin
`univariate` helper and flag it upstream.** A batch-invert over a packed
coefficient buffer is arguably `fgf`'s object (it is a field-vector kernel), and
`fgf`'s `ops` is where the AXPY/scatter kernels live. But `fgf` does not ship a
batch-invert today, and blocking the crate on an upstream addition is worse than
a small, clearly-flagged helper that migrates when `fgf` grows one.

*Revisit when* `fgf` exposes a `batch_invert` (or `ops::mul_inv_batch`); the
`univariate` helper is deleted and the call sites move upstream, the same way the
Cauchy inverse migrated from `srs` to `gfm`.

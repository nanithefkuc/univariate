# Rename fallout: `poly` → `univariate`

Noted at planning time (2026-08-16). The dependency-graph node currently labelled
`poly` is renamed `univariate` because the crate name `poly` is already taken
(the reservation under that name is not usable as a published crate). This file
is the change list; the root `README.md` and `AGENTS.md` already apply it.

## 1. Why the rename

The root `README.md` dependency graph and the `AGENTS.md` crate index both
carried a node/row named `poly` for "Univariate polynomial arithmetic over
GF(2^m)". That name is not available as a published crate. The permanent name
for the object — the univariate polynomial ring over GF(2^m) — is `univariate`.
The crate's `.plans` set is written directly under `univariate/` and uses the
new name throughout; this rename doc exists for parity with `gfm/.plans/09-fgf-rename.md`
and so the graph change is auditable in one place.

## 2. Root `README.md` — both notations

Applied (2026-08-16) at the rename; recorded here so a future revert or re-check
has the diff.

**Mermaid (`graph TD` block):**

- Node `poly["poly<br/><i>(Univariate Polynomial Arithmetic)</i>"]` →
  `univariate["univariate<br/><i>(Univariate Polynomial Arithmetic)</i>"]`.
- Edges, source side only:
  - `simdispatch --> poly` → `simdispatch --> univariate`
  - `fgf --> poly` → `fgf --> univariate`
  - `poly --> hasse` → `univariate --> hasse`
  - `poly --> funcfield` → `univariate --> funcfield`
  - `poly --> syndrome_engine` → `univariate --> syndrome_engine`
  - `poly --> gs_engine` → `univariate --> gs_engine`

**D2 (`direction: down` block):**

- Node `poly: "poly\n(Univariate Polynomial Arithmetic)" { class: impl }` →
  `univariate: "univariate\n(Univariate Polynomial Arithmetic)" { class: impl }`.
- Edges:
  - `L0.simdispatch -> L1.poly` → `L0.simdispatch -> L1.univariate`
  - `L1.fgf -> L1.poly` → `L1.fgf -> L1.univariate`
  - `L1.poly -> L1.hasse` → `L1.univariate -> L1.hasse`
  - `L1.poly -> L1.funcfield` → `L1.univariate -> L1.funcfield`
  - `L1.poly -> L2.syndrome_engine` → `L1.univariate -> L2.syndrome_engine`
  - `L1.poly -> L2.gs_engine` → `L1.univariate -> L2.gs_engine`

No other edges touched. The Wave-1 prose list in `README.md` does not name
`poly` explicitly (it enumerates `simdispatch`, `fgf`, `butterfly-fft`,
`sgraph`, `lattica`, `gfm`), so the prose needs no change; the node lives only
in the graph.

## 3. Root `AGENTS.md`

- The Layer code-fence list (the `L1 … mathematical domains` line) — `poly` →
  `univariate`.
- The crate-index table row — `| poly | 1 | 1 | Univariate polynomial
  arithmetic over GF(2^m) | planned | — |` → `| univariate (ex poly) | 1 | 1 |
  Univariate polynomial arithmetic over GF(2^m) | planned | — |`. The `(ex poly)`
  annotation matches the `gfm (ex linear-alg)` / `systematic-rs (ex srs)`
  precedent and records the provenance; it is removed when the `(ex …)` clutter
  is cleaned across the index, not before.

## 4. What is *not* renamed

- `gfm/src/poly.rs` — that is `gfm`'s *weak-Popov polynomial-row reduction*
  module, a different object (matrices of polynomials, linear algebra). It is
  not the univariate-ring node and does not move. The name collision is exactly
  why the L1 node cannot be `poly`: two things in the tree already answer to it.
- `gfm/.plans/02-architecture.md`'s planned `poly/` submodule (hankel /
  weak_popov / matrix) — `gfm`'s internal module namespace, unrelated to the L1
  crate name. See `08-risks.md` D1 for the seam this creates.
- Every `PolynomialError`, `Polynomial<F>`, `poly/` *module path inside
  `gs-engine`* — those are type/module names, not crate names; they migrate when
  `gs-engine` becomes a consumer (Phase 6), at which point they refer to
  `univariate::Polynomial`.

## 5. Smaller items

- This `.plans` set is written under `univariate/.plans/`; there is no
  `poly/.plans/` directory to migrate. Greenfield, like `gfm`'s was.
- The crate's own `README.md` (ground rule 4, written at scaffold time) opens
  with the AI-authorship `> [!WARNING]` header, not with the rename rationale.
  Provenance belongs in `CHANGELOG.md` and PRs, per the root style guide.
- `simdispatch --> univariate` (the L0→L1 edge every L1 node carries) is kept;
  `univariate` consumes the `Backend` re-export through `fgf`, not a direct
  `simdispatch` dependency — same as every other L1 node.

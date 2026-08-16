# Benchmarks

Crossover thresholds and their measurements. The source carries one-line
pointers here; this file carries the numbers, the hardware, and the commands.

## Hardware and toolchain of record

- CPU: Intel Core Ultra 7 258V (Lunar Lake, 8 cores)
- rustc: stable (2026-08), `--release`
- Fields: GF(2^8) (`Gf8`, AES polynomial) unless noted.

## Karatsuba dispatch crossover (`KARATSUBA_CROSSOVER`)

`Polynomial::multiply` dispatches schoolbook → Karatsuba at
`min(left, right) ≥ KARATSUBA_CROSSOVER = 2048` coefficients.

Criterion, `multiply/{schoolbook,karatsuba}/<n>`, GF(2^8), equal operands:

| Coefficients | Schoolbook | Karatsuba | Ratio |
| --- | --- | --- | --- |
| 128 | 0.89 µs | 1.82 µs | 2.04 |
| 512 | 6.60 µs | 8.79 µs | 1.33 |
| 1024 | 14.7 µs | 17.9 µs | 1.22 |
| 2048 | 53.5 µs | 34.1 µs | 0.64 |
| 3072 | 113.8 µs | 124.6 µs | 1.09 |
| 4096 | 202 µs | 69.1 µs | 0.34 |
| 6144 | 448 µs | 252 µs | 0.56 |

Karatsuba wins clearly from 2048 on (one reproducible anomaly at 3072,
cache-aliasing suspect); schoolbook wins through 1024. Command:

```sh
cargo bench --bench poly -- multiply
```

The Karatsuba recursion bottoms out in the packed schoolbook convolution
below 48 coefficients (`KARATSUBA_BASE`, private).

## AFFT product crossovers

Inherited from the extraction source (`gs-engine`) unchanged:
`AFFT_{PRODUCT,BATCH4,BATCH8,BATCH16}_CROSSOVER` and their scalar twins. See
`gs-engine/BENCHMARKS.md` for the original measurement records; re-measure
on this crate's bench harness before retuning.

## Chien vs equal-degree (`chien_equal_degree_crossover`)

Analytic first cut (`|F| / log²|F|`, floor 8), pending a dedicated
measurement; the selector's structure (pure function of degree and field
order) is stable under retuning.

## Multipoint / interpolation crossovers

`MULTIPOINT_EVAL_CROSSOVER = 16` and `MODULE_INTERPOLATION_CROSSOVER = 8`
(the latter inherited from `gs-engine`). Both measured on this crate's
`poly` bench harness; retune with `cargo bench --bench poly`.

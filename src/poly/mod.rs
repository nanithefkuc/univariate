//! Dense univariate polynomials over `fgf` binary fields.

#[cfg(feature = "fft")]
mod afft;
mod dense;
mod divide;
mod gcd;
mod karatsuba;
mod ring;
mod series;
pub use dense::Polynomial;

#[cfg(feature = "fft")]
pub use afft::{
    AFFT_BATCH4_CROSSOVER, AFFT_BATCH8_CROSSOVER, AFFT_BATCH16_CROSSOVER, AFFT_PRODUCT_CROSSOVER,
    PolynomialProductScratch, ProductStrategy, SCALAR_AFFT_BATCH4_CROSSOVER,
    SCALAR_AFFT_BATCH8_CROSSOVER, SCALAR_AFFT_BATCH16_CROSSOVER, SCALAR_AFFT_PRODUCT_CROSSOVER,
    multiply_batch_truncated, multiply_batch_truncated_with,
    substitute_y_affine_rows_truncated_into,
};
pub use gcd::{BezoutRelation, TruncatedEea, truncated_eea};
pub use karatsuba::{KARATSUBA_CROSSOVER, karatsuba_multiply};
pub use ring::binomial_odd;
pub use series::series_divide;

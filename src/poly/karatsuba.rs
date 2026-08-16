//! The Karatsuba middle tier of the product ladder.
//!
//! Between the schoolbook base and the AFFT batched product, operands are too
//! large for the quadratic convolution but too small to amortize a transform.
//! Karatsuba splits each operand at the half-degree, forms three half-size
//! products, and recombines with two additions; in characteristic two every
//! subtraction is an addition, so the classical sign bookkeeping vanishes.
//!
//! The recursion operates directly on packed coefficient byte slices with
//! depth-indexed scratch slots, so the tier costs field work and no
//! per-level allocation.

use alloc::vec::Vec;

use fgf::field::Elem as _;
use fgf::kernel::FieldKernels;
use fgf::ops;

use crate::error::{ConfigError, PolynomialError};

use super::dense::Polynomial;

/// Operand coefficient count at or above which [`Polynomial::multiply`]
/// dispatches to Karatsuba instead of schoolbook. Measured on GF(2^8)
/// (`BENCHMARKS.md`): schoolbook wins through 1024 coefficients (1.2-2x),
/// Karatsuba wins from 2048 (0.64x at 2048, 0.34x at 4096, 0.56x at 6144,
/// with a marginal 1.09x anomaly at 3072).
pub const KARATSUBA_CROSSOVER: usize = 2048;

/// Recursion base: below this shorter-operand size the recursion falls to
/// the packed schoolbook convolution. Measured; see `BENCHMARKS.md`.
const KARATSUBA_BASE: usize = 48;

/// One depth level of Karatsuba scratch: two split-sum operands and three
/// product slots, all byte buffers.
#[derive(Debug)]
struct KaratsubaScratch {
    sums: Vec<Vec<u8>>,
    products: Vec<Vec<u8>>,
}

impl KaratsubaScratch {
    const fn new() -> Self {
        Self {
            sums: Vec::new(),
            products: Vec::new(),
        }
    }

    /// Take the `index`-th sum slot, sized to at least `bytes` and zeroed
    /// over the active region.
    fn take_sum(&mut self, index: usize, bytes: usize) -> Vec<u8> {
        while self.sums.len() <= index {
            self.sums.push(Vec::new());
        }
        let mut slot = core::mem::take(&mut self.sums[index]);
        if slot.len() < bytes {
            slot.resize(bytes, 0);
        }
        slot[..bytes].fill(0);
        slot
    }

    /// Take the `index`-th product slot, sized to exactly `bytes` and zeroed.
    fn take_product(&mut self, index: usize, bytes: usize) -> Vec<u8> {
        while self.products.len() <= index {
            self.products.push(Vec::new());
        }
        let mut slot = core::mem::take(&mut self.products[index]);
        if slot.len() < bytes {
            slot.resize(bytes, 0);
        } else {
            slot.truncate(bytes);
        }
        slot.fill(0);
        slot
    }
}

/// Return the Karatsuba product of two nonzero polynomials.
///
/// Byte-identical to the schoolbook product; only the cost differs. The
/// recursion bottoms out in the packed schoolbook convolution once the
/// shorter operand falls below its measured base size.
///
/// # Errors
///
/// Returns [`PolynomialError::Config`] when the product buffer cannot be
/// reserved.
///
/// # Panics
///
/// The coefficient-aligned product expectation holds for the byte-exact
/// recombination this function computes itself.
pub fn karatsuba_multiply<F: FieldKernels>(
    a: &Polynomial<F>,
    b: &Polynomial<F>,
) -> Result<Polynomial<F>, PolynomialError> {
    let left = a.coefficient_count();
    let right = b.coefficient_count();
    debug_assert!(left > 0 && right > 0);
    let output_count = left
        .checked_add(right)
        .and_then(|sum| sum.checked_sub(1))
        .ok_or(ConfigError::GeometryOverflow {
            context: "polynomial product coefficients",
        })?;
    let mut destination = Vec::new();
    destination
        .try_reserve_exact(output_count * F::BYTES)
        .map_err(|_| ConfigError::AllocationFailed {
            context: "Karatsuba product",
            elements: output_count,
            element_size: F::BYTES,
        })?;
    destination.resize(output_count * F::BYTES, 0);
    let mut scratch = KaratsubaScratch::new();
    karatsuba_into::<F>(
        &mut destination,
        a.as_packed(),
        b.as_packed(),
        0,
        &mut scratch,
    );
    Ok(Polynomial::from_packed(destination).expect("coefficient-aligned product"))
}

/// Accumulate the packed schoolbook convolution `a · b` into zeroed `dst`.
fn schoolbook_into<F: FieldKernels>(dst: &mut [u8], a: &[u8], b: &[u8]) {
    debug_assert_eq!(dst.len() + F::BYTES, a.len() + b.len());
    for (index, chunk) in b.chunks_exact(F::BYTES).enumerate() {
        let scale = F::read(chunk);
        if scale.is_zero() {
            continue;
        }
        let offset = index * F::BYTES;
        let target = &mut dst[offset..offset + a.len()];
        ops::mul_add::<F>(target, scale, a);
    }
}

/// Byte XOR of `source` into `destination` (field addition in GF(2^m)).
fn xor_into(destination: &mut [u8], offset: usize, source: &[u8]) {
    for (byte, input) in destination[offset..].iter_mut().zip(source) {
        *byte ^= input;
    }
}

/// Write the Karatsuba product of packed `a` and `b` into zeroed `dst`.
fn karatsuba_into<F: FieldKernels>(
    dst: &mut [u8],
    a: &[u8],
    b: &[u8],
    depth: usize,
    scratch: &mut KaratsubaScratch,
) {
    let a_count = a.len() / F::BYTES;
    let b_count = b.len() / F::BYTES;
    if a_count.min(b_count) < KARATSUBA_BASE {
        schoolbook_into::<F>(dst, a, b);
        return;
    }

    let split = a_count.max(b_count).div_ceil(2);
    let split_bytes = split * F::BYTES;
    let a_split = split_bytes.min(a.len());
    let b_split = split_bytes.min(b.len());
    let (a_low, a_high) = a.split_at(a_split);
    let (b_low, b_high) = b.split_at(b_split);

    // z0 and z2 recurse on the external operand slices.
    let low_len = a_low.len() + b_low.len().saturating_sub(F::BYTES);
    let high_len = a_high.len() + b_high.len().saturating_sub(F::BYTES);
    let mut z0 = scratch.take_product(3 * depth, low_len);
    karatsuba_into::<F>(&mut z0, a_low, b_low, depth + 1, scratch);
    let mut z2 = scratch.take_product(3 * depth + 1, high_len);
    karatsuba_into::<F>(&mut z2, a_high, b_high, depth + 1, scratch);

    // Split sums, zero-padded to `split` coefficients each.
    let mut sum_a = scratch.take_sum(2 * depth, split_bytes);
    sum_a[..a_low.len()].copy_from_slice(a_low);
    xor_elementwise::<F>(&mut sum_a, a_high);
    let mut sum_b = scratch.take_sum(2 * depth + 1, split_bytes);
    sum_b[..b_low.len()].copy_from_slice(b_low);
    xor_elementwise::<F>(&mut sum_b, b_high);

    // middle = (a_low + a_high)(b_low + b_high), or empty when either sum
    // vanished (a pure power-of-two split with an empty high part).
    let middle_len = if sum_a.len() + sum_b.len() >= F::BYTES {
        sum_a.len() + sum_b.len() - F::BYTES
    } else {
        0
    };
    let mut middle = scratch.take_product(3 * depth + 2, middle_len);
    let sums_empty = sum_a.iter().all(|byte| *byte == 0) || sum_b.iter().all(|byte| *byte == 0);
    if !sums_empty {
        karatsuba_into::<F>(&mut middle, &sum_a, &sum_b, depth + 1, scratch);
    }

    // result = z0 ⊕ X^s·(middle ⊕ z0 ⊕ z2) ⊕ X^{2s}·z2 — five per-byte XOR
    // passes, all from owned locals so no pass aliases `dst`. Out-of-range
    // tails are truncated per byte, which never affects in-range results.
    xor_into(dst, 0, &z0);
    xor_into(dst, split_bytes, &middle);
    xor_into(dst, split_bytes, &z0);
    xor_into(dst, split_bytes, &z2);
    xor_into(dst, 2 * split_bytes, &z2);

    scratch.products[3 * depth] = z0;
    scratch.products[3 * depth + 1] = z2;
    scratch.products[3 * depth + 2] = middle;
    scratch.sums[2 * depth] = sum_a;
    scratch.sums[2 * depth + 1] = sum_b;
}

/// XOR one packed operand into another (in-place field addition), matching
/// element boundaries from the start.
fn xor_elementwise<F: FieldKernels>(destination: &mut [u8], source: &[u8]) {
    for (left, right) in destination
        .chunks_exact_mut(F::BYTES)
        .zip(source.chunks_exact(F::BYTES))
    {
        for (byte, input) in left.iter_mut().zip(right) {
            *byte ^= input;
        }
    }
}

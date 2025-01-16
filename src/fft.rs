//! FFT stuff
//! provide conversions between `crate::poly::Poly<F>` and `crate::poly::FFTPoly`

// TODO: Consider NTTs instead, or smth? (swanky_field_fft?)

use crate::{
    modular::lift_centered,
    poly::{Complex, FFTPoly, Poly},
};

use swanky_field::PrimeFiniteField;
use tfhe_fft::fft128::{f128, Plan};

// TODO? Should we tailor the FFT used depending on the (size of) the modulus?
// using < 26 bits should enable 64-bit FFT without precision loss (52 significant digits)
// And going about 52 (?) might get some benefit from doing an NTT instead

/// The FFT engine, enabling transformations between Poly and FFTPoly,
/// multiplication of the FFTPoly corresponds to a negacyclic convolution
#[derive(Debug)]
pub(crate) struct FFTPlan {
    inner: Plan,
}

/// Convert a field element into an f128; assuming the element fits in a single limb
fn field_to_f128<F: PrimeFiniteField>(x: F) -> f128 {
    let limb = lift_centered(x);

    // The most significant part
    let f = limb as f64;
    // The least significant part, whatever we might have lost in conversion to f64
    // NOTE: not constant time

    // approx + e = limb
    let approx = f as i64;
    let e = if approx >= limb {
        -((approx - limb) as f64)
    } else {
        (limb - approx) as f64
    };

    f128(f, e)
}

/// Convert an f64 into a field element, used as part of the f128 conversion
fn f64_to_field<F: PrimeFiniteField>(x: f64, modulus: u64) -> F {
    // This gets kinda annoying with negatives potentially not fitting in an i64, but only u64
    // Hence we first do the rounding in f64s, which should be fine (since close to 0 rounds clearly)
    // We also once again introduce assumptions on limb size by using u128
    F::try_from_int::<1>(
        if x >= 0f64 {
            (x as u128 % modulus as u128) as u64
        } else {
            // Needs an extra mod when x0 = modulus
            // TODO: improve?
            (modulus - ((-x) as u128 % modulus as u128) as u64) % modulus
        }
        .into(),
    )
    .expect("Modular reduction gone wrong in f64_to_field")
}

/// Convert an f128 back into a field element; assuming the field modulus fits in a single limb
fn f128_to_field<F: PrimeFiniteField>(x: f128) -> F {
    // not ct, and ugly, but good enough

    let modulus = F::modulus_int::<1>().as_words()[0];
    let x0 = x.0.round();
    let x1 = x.1.round();
    f64_to_field::<F>(x0, modulus) + f64_to_field(x1, modulus)
}

impl FFTPlan {
    /// Create a new FFT plan taking degree `dim - 1` real (finite field) valued polynomials
    /// `dim` should be a power of two
    pub(crate) fn new(dim: usize) -> Self {
        debug_assert_eq!(dim.count_ones(), 1);

        Self {
            inner: Plan::new(dim >> 1),
        }
    }

    // TODO? Cached monomial transforms

    /// The forward FFT transformation
    pub(crate) fn fwd<F: PrimeFiniteField>(&self, poly: Poly<F>) -> FFTPoly {
        debug_assert_eq!(F::MIN_LIMBS_NEEDED, 1);
        debug_assert_eq!(poly.degree(), 2 * self.inner.fft_size() - 1);

        // I'm not *exactly* sure about the precision we get here
        // as it's going through f128
        // tfhe-rs uses a 64-bit modulus,
        // so we can probably get enough without too much noise growth

        let n = self.inner.fft_size();
        let mut re0 = vec![0f64; n];
        let mut re1 = vec![0f64; n];
        let mut im0 = vec![0f64; n];
        let mut im1 = vec![0f64; n];

        for i in 0..n {
            let a = field_to_f128(poly.get(i));
            let b = field_to_f128(poly.get(i + n));
            re0[i] = a.0;
            re1[i] = a.1;
            im0[i] = b.0;
            im1[i] = b.1;
        }

        self.inner.fwd(&mut re0, &mut re1, &mut im0, &mut im1);

        let mut res = Vec::with_capacity(n);
        for i in 0..n {
            res.push(Complex::new(f128(re0[i], re1[i]), f128(im0[i], im1[i])));
        }

        FFTPoly::new(res)
    }

    /// The inverse FFT transformation
    pub(crate) fn inv<F: PrimeFiniteField>(&self, poly: FFTPoly) -> Poly<F> {
        debug_assert_eq!(F::MIN_LIMBS_NEEDED, 1);
        debug_assert_eq!(poly.0.len(), self.inner.fft_size());

        // I'm not *exactly* sure about the precision we get here
        // as it's going through f128
        // tfhe-rs uses a 64-bit modulus,
        // so we can probably get enough without too much noise growth

        let n = self.inner.fft_size();
        let mut re0 = vec![0f64; n];
        let mut re1 = vec![0f64; n];
        let mut im0 = vec![0f64; n];
        let mut im1 = vec![0f64; n];

        let scale = 1f64 / (n as f64);
        for i in 0..n {
            re0[i] = poly.0[i].re.0 * scale;
            re1[i] = poly.0[i].re.1 * scale;
            im0[i] = poly.0[i].im.0 * scale;
            im1[i] = poly.0[i].im.1 * scale;
        }

        self.inner.inv(&mut re0, &mut re1, &mut im0, &mut im1);

        let mut res = Vec::with_capacity(2 * n);
        res.resize(2 * n, F::ZERO);
        for i in 0..n {
            res[i] = f128_to_field(f128(re0[i], re1[i]));
            res[i + n] = f128_to_field(f128(im0[i], im1[i]));
        }

        Poly(res)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::params::{Params, TESTTYPE};
    use crate::test_utils::*;

    use swanky_field::FiniteRing;

    // NOTE: All these tests currently rely on the modulus being small enough that
    // no precision is lost (at all, after rounding)
    const _SMALL_MOD: () =
        assert!(<<TESTTYPE as swanky_field::FiniteField>::NumberOfBitsInBitDecomposition as generic_array::typenum::Unsigned>::U8 <= 32);

    #[test]
    fn field_f128_roundtrip() {
        let rng = &mut rng(None);
        for _ in 0..100 {
            let x = TESTTYPE::random(rng);
            assert_eq!(x, f128_to_field(field_to_f128(x)));
        }
    }

    #[test]
    fn field_f128_mul_roundtrip() {
        let rng = &mut rng(None);
        for _ in 0..100 {
            let x = TESTTYPE::random(rng);
            let y = TESTTYPE::random(rng);
            let z: TESTTYPE = f128_to_field(field_to_f128(x) * field_to_f128(y));
            assert_eq!(x * y, z);
        }
    }

    #[test]
    fn field_f128_sub_roundtrip() {
        let rng = &mut rng(None);
        for _ in 0..100 {
            let x = TESTTYPE::random(rng);
            let y = TESTTYPE::random(rng);
            let z: TESTTYPE = f128_to_field(field_to_f128(x) - field_to_f128(y));
            assert_eq!(x - y, z);
        }
    }

    #[test]
    fn field_f128_mul_sub_roundtrip() {
        let rng = &mut rng(None);
        for _ in 0..100 {
            let x = TESTTYPE::random(rng);
            let y = TESTTYPE::random(rng);
            let a = TESTTYPE::random(rng);
            let b = TESTTYPE::random(rng);
            let z: TESTTYPE = f128_to_field(
                field_to_f128(a) * field_to_f128(b) - field_to_f128(x) * field_to_f128(y),
            );
            assert_eq!(a * b - x * y, z);
        }
    }

    #[test]
    fn f128_to_field_negative_modulus() {
        let modulus = TESTTYPE::modulus_int::<1>().as_words()[0];
        let f = f128(-(modulus as f64), 0.0);
        assert_eq!(f128_to_field::<TESTTYPE>(f), TESTTYPE::ZERO);
        let f = f128(0.0, -(modulus as f64));
        assert_eq!(f128_to_field::<TESTTYPE>(f), TESTTYPE::ZERO);
    }

    fn fft_test_roundtrip_size(n: usize) {
        let rng = &mut rng(None);
        // Needs to be >= 64 for tfhe_fft to be happy :/
        let plan = FFTPlan::new(n);
        let mut poly = Poly(vec![TESTTYPE::ZERO; n]);
        poly.iter_mut().for_each(|x| *x = TESTTYPE::random(rng));
        let fftpoly = plan.fwd(poly.clone());
        let result: Poly<TESTTYPE> = plan.inv(fftpoly);
        assert_eq!(poly, result);
    }

    #[test]
    fn fft_roundtrip_64() {
        fft_test_roundtrip_size(64);
    }

    #[test]
    fn fft_roundtrip() {
        fft_test_roundtrip_size(1 << crate::params::TestParams::LOG_DEG_NTRU);
    }

    fn fft_test_mult<F: PrimeFiniteField>(n: usize) {
        let rng = &mut rng(None);
        let plan = FFTPlan::new(n);

        for _ in 0..10 {
            let a = Poly::<F>::rand(n, rng);
            let b = Poly::<F>::rand(n, rng);

            let c_expected = Poly(slow_negacyclic_mult(&a.0, &b.0, F::ZERO));

            let a_fft = plan.fwd(a);
            let b_fft = plan.fwd(b);
            let c_fft = FFTPoly::new(Vec::from_iter(
                a_fft
                    .0
                    .into_iter()
                    .zip(b_fft.0.into_iter())
                    .map(|(x, y)| x * y),
            ));

            let c: Poly<F> = plan.inv(c_fft);
            assert!(c
                .0
                .into_iter()
                .zip(c_expected.0.into_iter())
                .all(|(x, y)| x == y))
        }
    }

    #[test]
    fn fft_mult_64() {
        fft_test_mult::<TESTTYPE>(64);
    }

    #[test]
    fn fft_mult() {
        fft_test_mult::<TESTTYPE>(1 << crate::params::TestParams::LOG_DEG_NTRU);
    }
}

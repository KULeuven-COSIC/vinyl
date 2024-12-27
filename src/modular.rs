//! Modular integer types and associated traits

use crate::params::Rng;
/// A type that can be modulo-switched onto `Target`
/// Functionally similar to [`Into`], but semantically more restricted
pub trait ModSwitch<Target> {
    fn modswitch(self) -> Target;
}

const _REQUIRE_64BIT_LIMB: () = assert!(crypto_bigint::Word::BITS == 64);

/// Integers modulo the 61-bit prime 2^61 - 1
pub type Q61 = swanky_field_f61p::F61p;

mod m31 {
    #![allow(clippy::suspicious_op_assign_impl)]

    use std::ops::{AddAssign, MulAssign, SubAssign};

    use crypto_bigint::subtle::{ConditionallySelectable, ConstantTimeEq};
    use swanky_field::{FiniteField, FiniteRing, PrimeFiniteField};
    use swanky_serialization::{BiggerThanModulus, CanonicalSerialize};

    // Mostly referenced from the F61p implementation
    /// Integers modulo the Mersenne prime 2^31 - 1
    #[allow(clippy::derived_hash_with_manual_eq)]
    #[derive(Copy, Clone, Debug, Eq, PartialOrd, Ord, Hash)]
    pub struct M31(u32);
    const MOD: u32 = (1 << 31) - 1;

    impl TryFrom<u128> for M31 {
        type Error = BiggerThanModulus;

        fn try_from(value: u128) -> Result<Self, Self::Error> {
            if value < MOD as u128 {
                Ok(Self(value as u32))
            } else {
                Err(BiggerThanModulus)
            }
        }
    }

    impl ConstantTimeEq for M31 {
        fn ct_eq(&self, other: &Self) -> crypto_bigint::subtle::Choice {
            self.0.ct_eq(&other.0)
        }
    }

    impl ConditionallySelectable for M31 {
        fn conditional_select(a: &Self, b: &Self, choice: crypto_bigint::subtle::Choice) -> Self {
            Self(u32::conditional_select(&a.0, &b.0, choice))
        }
    }

    fn reduce(x: u64) -> M31 {
        let y = (x & MOD as u64) + (x >> 31);
        M31(u32::conditional_select(
            &(y as u32),
            &(y.wrapping_sub(MOD as u64) as u32),
            ((y >= MOD as u64) as u8).into(),
        ))
    }

    impl AddAssign<&M31> for M31 {
        fn add_assign(&mut self, rhs: &Self) {
            *self = reduce(self.0 as u64 + rhs.0 as u64);
        }
    }

    impl SubAssign<&M31> for M31 {
        fn sub_assign(&mut self, rhs: &Self) {
            *self = reduce(self.0 as u64 + MOD as u64 - rhs.0 as u64);
        }
    }

    impl MulAssign<&M31> for M31 {
        fn mul_assign(&mut self, rhs: &Self) {
            *self = reduce(self.0 as u64 * rhs.0 as u64);
        }
    }

    impl CanonicalSerialize for M31 {
        type Serializer = swanky_serialization::ByteElementSerializer<Self>;
        type Deserializer = swanky_serialization::ByteElementDeserializer<Self>;
        type ByteReprLen = generic_array::typenum::U4;
        type FromBytesError = BiggerThanModulus;

        fn from_bytes(
            bytes: &generic_array::GenericArray<u8, Self::ByteReprLen>,
        ) -> Result<Self, Self::FromBytesError> {
            let buf = <[u8; 4]>::from(*bytes);
            let x = u32::from_le_bytes(buf);
            Self::try_from(x as u128)
        }

        fn to_bytes(&self) -> generic_array::GenericArray<u8, Self::ByteReprLen> {
            self.0.to_le_bytes().into()
        }
    }

    swanky_field::field_ops!(M31);

    impl FiniteRing for M31 {
        fn from_uniform_bytes(x: &[u8; 16]) -> Self {
            Self((u128::from_le_bytes(*x) % MOD as u128) as u32)
        }

        fn random<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
            Self(rng.gen_range(0..MOD))
        }

        const ZERO: Self = Self(0);
        const ONE: Self = Self(1);
    }

    impl FiniteField for M31 {
        type PrimeField = Self;

        const GENERATOR: Self = Self(7);

        type NumberOfBitsInBitDecomposition = generic_array::typenum::U31;

        fn bit_decomposition(
            &self,
        ) -> generic_array::GenericArray<bool, Self::NumberOfBitsInBitDecomposition> {
            swanky_field::standard_bit_decomposition(self.0 as u128)
        }

        fn inverse(&self) -> Self {
            if *self == Self::ZERO {
                panic!("Zero cannot be inverted");
            }
            self.pow_var_time(MOD as u128 - 2)
        }

        fn polynomial_modulus() -> swanky_field::polynomial::Polynomial<Self::PrimeField> {
            swanky_field::polynomial::Polynomial::x()
        }
    }

    impl PrimeFiniteField for M31 {
        fn modulus_int<const LIMBS: usize>() -> crypto_bigint::Uint<LIMBS> {
            crypto_bigint::Uint::<LIMBS>::from_u32(MOD)
        }

        fn into_int<const LIMBS: usize>(&self) -> crypto_bigint::Uint<LIMBS> {
            crypto_bigint::Uint::<LIMBS>::from_u32(self.0)
        }

        fn try_from_int<const LIMBS: usize>(
            x: crypto_bigint::Uint<LIMBS>,
        ) -> crypto_bigint::subtle::CtOption<Self> {
            use crypto_bigint::subtle::ConstantTimeLess;

            crypto_bigint::subtle::CtOption::new(
                Self(x.as_words()[0] as u32),
                x.ct_lt(&Self::modulus_int()),
            )
        }
    }
}

pub use m31::M31;

/// Sample an element from the discrete Gaussian distribution
/// centered around zero, with given standard deviation.
pub(crate) fn sample_discrete_gaussian(stdev: f64, rng: &mut impl Rng) -> u64 {
    // TODO? Do this properly instead of through rounding
    let distr = rand_distr::Normal::new(0f64, stdev)
        .expect("Invalid distribution to sample discrete Gaussian");
    rng.sample(distr).round() as u64
}

/// Lift a bit into the finite field. The bit must be 0 or 1.
pub(crate) fn bit_to_field<F: swanky_field::FiniteField>(bit: u8) -> F {
    debug_assert!(bit == 0 || bit == 1);

    F::conditional_select(&F::ZERO, &F::ONE, bit.into())
}

/// Lift an int into the field.
/// This only works for "small" integers, as we assume it fits in a single crypto_bigint limb
///
/// # Panics
///
/// Panics when the value is not properly reduced to [0, modulus)
// TODO: Make this general enough for larger Uints?
pub(crate) fn int_to_field<F: swanky_field::PrimeFiniteField>(val: crypto_bigint::Limb) -> F {
    F::try_from_int::<1>(val.into()).expect("Value not reduce for field modulus")
}

pub(crate) fn lift_centered<F: swanky_field::PrimeFiniteField>(val: F) -> i64 {
    let modulus: i64 = F::modulus_int::<1>().as_words()[0]
        .try_into()
        .expect("Needs modulus to fit in i64");

    let lift = val.into_int::<1>().as_words()[0].try_into().unwrap();
    if lift < modulus / 2 {
        lift
    } else {
        lift - modulus
    }
}

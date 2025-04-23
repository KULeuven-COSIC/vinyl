//! Modular integer types and associated traits

use crate::params::Rng;

/// A type that can be modulo-switched onto `Target`
/// Functionally similar to [`Into`], but semantically more restricted
pub trait ModSwitch<Target> {
    fn modswitch(self) -> Target;
}

fn randomized_rounding_div(a: u128, b: u128) -> u128 {
    use rand::Rng;

    let q = a / b;
    let r = a % b;
    let round = rand::thread_rng().gen_range(0..b);
    q + (round <= r) as u128
}

// TODO: the usual, assumes everything fits in a limb
impl<F1, F2> ModSwitch<F2> for F1
where
    // A lot of these bounds are because annoyance from ModReduce :/
    F1: Modular,
    F1::Single: Into<u128>,
    F2: Modular + ModReduce,
    F2::Double: From<F2::Single> + std::ops::Rem<Output = F2::Double>,
    F2::Single: TryFrom<F2::Double>,
    <F2::Single as TryFrom<F2::Double>>::Error: std::fmt::Debug,
    F2::Single: Into<u128>,
    F2::Double: TryFrom<u128>,
    <F2::Double as TryFrom<u128>>::Error: std::fmt::Debug,
{
    fn modswitch(self) -> F2 {
        let mod_from = F1::MOD.into();
        let mod_to = F2::MOD.into();
        let n = self.extract().into();
        F2::reduce(
            randomized_rounding_div(mod_to * n, mod_from)
                .try_into()
                .expect("u128 doesn't fit in Modular::Double after modswitch"),
        )
    }
}

const _REQUIRE_64BIT_LIMB: () = assert!(crypto_bigint::Word::BITS == 64);

/// A type that has a modulus
pub trait Modular {
    /// The type needed to represent a fully reduced number
    type Single;
    /// Type type needed to represent the product of two fully reduced numbers
    type Double;

    const MOD: Self::Single;

    /// Create a new object, assuming it's smaller than the modulus
    fn new_unchecked(x: Self::Single) -> Self;

    /// Get the underlying value
    fn extract(&self) -> Self::Single;
}

pub(crate) trait ModReduce
where
    Self: Modular + Sized,
    // Ideally, the next bounds would only be there for the default implementation of `reduce`
    Self::Double: From<Self::Single> + std::ops::Rem<Output = Self::Double>,
    Self::Single: TryFrom<Self::Double>,
    <Self::Single as TryFrom<Self::Double>>::Error: std::fmt::Debug,
{
    #[inline]
    fn reduce(val: Self::Double) -> Self {
        Self::new_unchecked(
            Self::Single::try_from(val % Self::Double::from(Self::MOD))
                .expect("Modular reduction should always fit in Single"),
        )
    }
}

pub(crate) mod smallfields {
    macro_rules! small_field {
        ($mod_name:ident::$name:ident, $single:ty, $double:ty, $mod:expr, $generator:literal) => {
            impl $crate::modular::ModReduce for $mod_name::$name { }
            $crate::modular::smallfields::small_field!(@custom_reduce, $mod_name::$name, $single, $double, $mod, $generator);
        };

        (@custom_reduce, $mod_name:ident::$name:ident, $single:ty, $double:ty, $mod:expr, $generator:literal) => {
            mod $mod_name {
                #![allow(clippy::suspicious_op_assign_impl)]
                use $crate::modular::{Modular, ModReduce};
                use swanky_field::{FiniteRing, FiniteField};
                use crypto_bigint::subtle::ConstantTimeEq;

                #[allow(clippy::derived_hash_with_manual_eq)]
                #[derive(Copy, Clone, Debug, Eq, PartialOrd, Ord, Hash)]
                pub struct $name($single);

                impl Modular for $name {
                    type Single = $single;
                    type Double = $double;

                    const MOD: Self::Single = $mod;

                    #[inline]
                    fn new_unchecked(x: Self::Single) -> Self {
                        Self(x)
                    }

                    #[inline]
                    fn extract(&self) -> Self::Single {
                        self.0
                    }
                }

                impl TryFrom<u128> for $name where Self: Modular {
                    type Error = swanky_serialization::BiggerThanModulus;

                    #[inline]
                    fn try_from(value: u128) -> Result<Self, Self::Error> {
                        if value < Self::MOD.into() {
                            Ok(Self(value.try_into().expect("Value smaller than modulus that fits")))
                        } else {
                            Err(swanky_serialization::BiggerThanModulus)
                        }
                    }
                }

                impl ConstantTimeEq for $name {
                    #[inline]
                    fn ct_eq(&self, other: &Self) -> crypto_bigint::subtle::Choice {
                        self.0.ct_eq(&other.0)
                    }
                }

                impl crypto_bigint::subtle::ConditionallySelectable for $name {
                    #[inline]
                    fn conditional_select(a: &Self, b: &Self, choice: crypto_bigint::subtle::Choice) -> Self {
                        Self(<$single>::conditional_select(&a.0, &b.0, choice))
                    }
                }

                impl std::ops::AddAssign<&$name> for $name {
                    #[inline]
                    fn add_assign(&mut self, rhs: &Self) {
                        *self = Self::reduce(self.0 as $double + rhs.0 as $double);
                    }
                }

                impl std::ops::SubAssign<&$name> for $name {
                    #[inline]
                    fn sub_assign(&mut self, rhs: &Self) {
                        *self = Self::reduce(self.0 as $double + Self::MOD as $double - rhs.0 as $double);
                    }
                }

                impl std::ops::MulAssign<&$name> for $name {
                    #[inline]
                    fn mul_assign(&mut self, rhs: &Self) {
                        *self = Self::reduce(self.0 as $double * rhs.0 as $double);
                    }
                }

                impl swanky_serialization::CanonicalSerialize for $name {
                    type Serializer = swanky_serialization::ByteElementSerializer<Self>;
                    type Deserializer = swanky_serialization::ByteElementDeserializer<Self>;
                    type ByteReprLen = <generic_array::typenum::Const<{ <$single>::BITS as usize / 8 }> as generic_array::IntoArrayLength>::ArrayLength;
                    type FromBytesError = swanky_serialization::BiggerThanModulus;

                    #[inline]
                    fn from_bytes(
                        bytes: &generic_array::GenericArray<u8, Self::ByteReprLen>,
                    ) -> Result<Self, Self::FromBytesError> {
                        let buf = <[u8; <$single>::BITS as usize / 8]>::from(*bytes);
                        let x = <$single>::from_le_bytes(buf);
                        Self::try_from(x as u128)
                    }

                    #[inline]
                    fn to_bytes(&self) -> generic_array::GenericArray<u8, Self::ByteReprLen> {
                        self.0.to_le_bytes().into()
                    }
                }

                impl swanky_field::FiniteRing for $name {
                    #[inline]
                    fn from_uniform_bytes(x: &[u8; 16]) -> Self {
                        Self((u128::from_le_bytes(*x) % Self::MOD as u128) as $single)
                    }

                    #[inline]
                    fn random<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                        Self(rng.gen_range(0..Self::MOD))
                    }

                    const ZERO: Self = Self(0);
                    const ONE: Self = Self(1);
                }

                impl swanky_field::FiniteField for $name {
                    type PrimeField = Self;

                    const GENERATOR: Self = Self($generator);

                    type NumberOfBitsInBitDecomposition = <generic_array::typenum::Const<{ <$single>::BITS as usize - <Self as Modular>::MOD.leading_zeros() as usize }> as generic_array::IntoArrayLength>::ArrayLength;

                    #[inline]
                    fn bit_decomposition(
                        &self,
                    ) -> generic_array::GenericArray<bool, Self::NumberOfBitsInBitDecomposition> {
                        swanky_field::standard_bit_decomposition(self.0 as u128)
                    }

                    // TODO? Faster
                    #[inline]
                    fn inverse(&self) -> Self {
                        if *self == Self::ZERO {
                            panic!("Zero cannot be inverted");
                        }
                        self.pow_var_time(Self::MOD as u128 - 2)
                    }

                    #[inline]
                    fn polynomial_modulus() -> swanky_field::polynomial::Polynomial<Self::PrimeField> {
                        swanky_field::polynomial::Polynomial::x()
                    }
                }

                impl swanky_field::PrimeFiniteField for $name {
                    #[inline]
                    fn modulus_int<const LIMBS: usize>() -> crypto_bigint::Uint<LIMBS> {
                        crypto_bigint::Uint::<LIMBS>::from(Self::MOD)
                    }

                    #[inline]
                    fn into_int<const LIMBS: usize>(&self) -> crypto_bigint::Uint<LIMBS> {
                        crypto_bigint::Uint::<LIMBS>::from(self.0)
                    }

                    #[inline]
                    fn try_from_int<const LIMBS: usize>(
                        x: crypto_bigint::Uint<LIMBS>,
                    ) -> crypto_bigint::subtle::CtOption<Self> {
                        use crypto_bigint::subtle::ConstantTimeLess;

                        crypto_bigint::subtle::CtOption::new(
                            // The condition should make sure the default unwrap doesn't matter
                            // Probably not ct though...
                            Self(x.as_words()[0].try_into().unwrap_or_default()),
                            x.ct_lt(&Self::modulus_int()),
                        )
                    }
                }

                swanky_field::field_ops!($name);
            }
            #[allow(unused_imports)]
            pub use $mod_name::$name;
        };
    }

    pub(crate) use small_field;
}

smallfields::small_field!(q20::Q20, u32, u64, 912829, 2);
smallfields::small_field!(q17::Q17, u32, u64, 92683, 3);
smallfields::small_field!(q30::Q30, u32, u64, 1073741827, 2);
smallfields::small_field!(q35::Q35, u64, u128, 34359738421, 2);

impl ModReduce for M31 {
    #[inline]
    fn reduce(val: Self::Double) -> Self {
        use crypto_bigint::subtle::ConditionallySelectable;

        let y = (val & Self::MOD as u64) + (val >> 31);
        Self::new_unchecked(u32::conditional_select(
            &(y as u32),
            &(y.wrapping_sub(Self::MOD as u64) as u32),
            ((y >= Self::MOD as u64) as u8).into(),
        ))
    }
}

smallfields::small_field!(@custom_reduce, m31::M31, u32, u64, (1 << 31) - 1, 7);

impl ModReduce for M61 {
    #[inline]
    fn reduce(val: Self::Double) -> Self {
        use crypto_bigint::subtle::ConditionallySelectable;

        let y = (val & Self::MOD as u128) + (val >> 61);
        Self::new_unchecked(u64::conditional_select(
            &(y as u64),
            &(y.wrapping_sub(Self::MOD as u128) as u64),
            ((y >= Self::MOD as u128) as u8).into(),
        ))
    }
}
smallfields::small_field!(@custom_reduce, m61::M61, u64, u128, (1 << 61) - 1, 37);

#[derive(Debug)]
pub struct Z2k<const K: usize>(u32);

impl<const K: usize> Z2k<K> {
    const _K_FITS: () = assert!(K < 32);
    pub const MASK: u32 = (1 << K) - 1;
}

impl<const K: usize> Modular for Z2k<K> {
    type Single = u32;
    type Double = u64;
    const MOD: u32 = 1 << K;

    fn new_unchecked(x: Self::Single) -> Self {
        Self(x)
    }

    fn extract(&self) -> u32 {
        self.0
    }
}

impl<const K: usize> ModReduce for Z2k<K> {
    fn reduce(val: Self::Double) -> Self {
        Self::new_unchecked((val & Self::MASK as u64) as u32)
    }
}

/// Sample an element from the discrete Gaussian distribution
/// centered around zero, with given standard deviation.
pub(crate) fn sample_discrete_gaussian(stdev: f64, rng: &mut impl Rng) -> u64 {
    // TODO? Do this properly instead of through rounding
    let distr = rand_distr::Normal::new(0f64, stdev)
        .expect("Invalid distribution to sample discrete Gaussian");
    rng.sample(distr).round() as u64
}

/// Sample a random ternary element, embedded into the field F
#[cfg_vis::cfg_vis(feature = "bench", pub)]
pub(crate) fn sample_ternary<F: swanky_field::PrimeFiniteField>(rng: &mut impl Rng) -> F {
    let r = rng.gen_range(0..3u8).into();
    int_to_field::<F>(r) - F::ONE
}

pub(crate) fn sample_gaussian_ternary<F: swanky_field::PrimeFiniteField>(rng: &mut impl Rng) -> F {
    int_to_field::<F>(rng.gen_range(0..=1u8).into())
        - int_to_field::<F>(rng.gen_range(0..=1u8).into())
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
    F::try_from_int::<1>(val.into()).expect("Value not reduced for field modulus")
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::{params::TESTTYPE, test_utils::*};

    use swanky_field::FiniteRing;

    #[test]
    fn modswitch_self() {
        let rng = &mut rng(None);

        for _ in 0..100 {
            let x = TESTTYPE::random(rng);
            assert_eq!(x, x.modswitch());
        }
    }

    #[test]
    fn modswitch_eq_modulus() {
        for _ in 0..100 {
            let x = M61::ZERO - M61::ONE;
            let y: M31 = x.clone().modswitch();
            // we care about the line above not panicking
            core::hint::black_box(y);
        }
    }

    #[test]
    fn modswitch_field() {
        let rng = &mut rng(None);

        for _ in 0..100 {
            let x = TESTTYPE::random(rng);
            //TODO
        }
    }
}

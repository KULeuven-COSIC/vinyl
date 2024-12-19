//! Modular integer types and associated traits

use crate::params::Rng;

/// A type that can be modulo-switched onto `Target`
/// Functionally similar to [`Into`], but semantically more restricted
pub trait ModSwitch<Target> {
    fn modswitch(self) -> Target;
}

/// Integers modulo the 64-bit prime 2^64 - 59
pub type Q61 = swanky_field_f61p::F61p;

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

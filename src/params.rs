//! Parameters defining the entire scheme

use crate::modular::Q61;

/// A cryptographically secure RNG trait for convenience
pub trait Rng: rand::Rng + rand::CryptoRng {}
impl<T: rand::Rng + rand::CryptoRng> Rng for T {}

/// Parameters that define the size, security and efficiency of the entire FINAL FHE scheme
#[derive(Debug)]
pub struct Params<BootInt, BaseInt> {
    _tp: std::marker::PhantomData<(BootInt, BaseInt)>,
    /// The LWE dimension, for binary LWE secrets
    pub(crate) dim_lwe: usize,
    /// log2(N) where the NTRU ciphertext are defined over mod X^N + 1
    pub(crate) log_deg_ntru: u64,
    /// The number of components for gadget decomposition
    pub(crate) dim_ngs: u64,
    /// The (discrete) distribution used to sample LWE noise
    pub(crate) err_stdev_lwe: f64,
    /// The (discrete) distribution used for fresh NTRU noise polynomials
    pub(crate) err_stdev_ntru: f64,
}

// TODO: check!
/// Some parameters to test the scheme with during development
pub const TESTPARAMS: Params<Q61, Q61> = Params {
    _tp: std::marker::PhantomData,
    dim_lwe: 930,
    log_deg_ntru: 14,
    dim_ngs: 3,
    err_stdev_lwe: 4.39,
    err_stdev_ntru: 4.39,
};

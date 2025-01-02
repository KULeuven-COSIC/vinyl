//! Parameters defining the entire scheme

use crate::modular::{M31, M61};

use std::sync::LazyLock;

/// A cryptographically secure RNG trait for convenience
pub trait Rng: rand::Rng + rand::CryptoRng {}
impl<T: rand::Rng + rand::CryptoRng> Rng for T {}

/// Parameters that define the size, security and efficiency of the entire FINAL FHE scheme
#[derive(Debug)]
pub struct ParamsInner<BootInt, BaseInt> {
    _tp: std::marker::PhantomData<(BootInt, BaseInt)>,
    /// The LWE dimension, for binary LWE secrets
    pub(crate) dim_lwe: usize,
    /// log2(N) where the NTRU ciphertext are defined over mod X^N + 1
    pub(crate) log_deg_ntru: u64,
    /// The number of components for gadget decomposition
    pub(crate) dim_ngs: usize,
    /// The base for gadget decomposition
    pub(crate) gadget_base: BootInt,
    /// The (discrete) distribution used to sample LWE noise
    pub(crate) err_stdev_lwe: f64,
    /// The (discrete) distribution used for fresh NTRU noise polynomials
    pub(crate) err_stdev_ntru: f64,

    // NOTE: Putting the lazylock on this doesn't entirely solve
    // the issue of caching initialization, as (in theory) you could
    // use multiple `Params` instances with the same ntru degree and
    // hence duplicate FFTPlan work/space.
    // It's fine in practice though :)
    /// The fft plan/engine to use
    pub(crate) fft: crate::fft::FFTPlan,
}

/// Parameters that define the size, security and efficiency of the entire FINAL FHE scheme
pub type Params<BootInt, BaseInt> = LazyLock<ParamsInner<BootInt, BaseInt>>;

#[cfg(test)]
pub(crate) type TESTTYPE = M31;

// TODO: check!
/// Some parameters to test the scheme with during development
pub static TESTPARAMS: Params<M31, M31> = Params::new(|| ParamsInner {
    _tp: std::marker::PhantomData,
    dim_lwe: 930,
    log_deg_ntru: 10,
    dim_ngs: 3,
    gadget_base: <M31 as crate::modular::Modular>::new_unchecked(16),
    err_stdev_lwe: 4.39,
    err_stdev_ntru: 4.39,

    fft: crate::fft::FFTPlan::new(1 << 10),
});

// TODO: check, this is entirely random garbage
pub static LARGEPARAMS: Params<M61, M61> = Params::new(|| ParamsInner {
    _tp: std::marker::PhantomData,
    dim_lwe: 2048,
    log_deg_ntru: 14,
    dim_ngs: 4,
    gadget_base: <M61 as crate::modular::Modular>::new_unchecked(16),
    err_stdev_lwe: 12.0,
    err_stdev_ntru: 12.0,

    fft: crate::fft::FFTPlan::new(1 << 14),
});

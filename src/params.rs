//! Parameters defining the entire scheme

use crate::modular::{Mod2PowN, Q64};

/// Parameters that define the size, security and efficiency of the entire FINAL FHE scheme
#[derive(Debug)]
pub struct Params<BootInt, BaseInt, NTRUExpInt> {
    _tp: std::marker::PhantomData<(BootInt, BaseInt, NTRUExpInt)>,
    /// The LWE dimension, for binary LWE secrets
    pub(crate) dim_lwe: u64,
    /// log2(N) where the NTRU ciphertext are defined over mod X^N + 1
    pub(crate) log_deg_ntru: u64,
    /// The number of components for gadget decomposition
    pub(crate) dim_ngs: u64,
    /// (Discrete) Gaussian standard deviation for fresh LWE noise
    pub(crate) err_stdev_lwe: f64,
    /// (Discrete) Gaussian standard deviation for fresh NTRU noise polynomials
    pub(crate) err_stdev_ntru: f64,
}

// TODO: check!
/// Some parameters to test the scheme with during development
pub const TESTPARAMS: Params<Q64, Q64, Mod2PowN<15>> = Params {
    _tp: std::marker::PhantomData,
    dim_lwe: 930,
    log_deg_ntru: 14,
    dim_ngs: 3,
    err_stdev_lwe: 4.39,
    err_stdev_ntru: 4.39,
};

//! Parameters defining the entire scheme

use crate::modular::{Z2k, M31, M61};

use std::sync::LazyLock;

/// A cryptographically secure RNG trait for convenience
pub trait Rng: rand::Rng + rand::CryptoRng {}
impl<T: rand::Rng + rand::CryptoRng> Rng for T {}

/// Parameters that define the size, security and efficiency of the entire FINAL FHE scheme
#[derive(Debug)]
pub struct ParamsInner<BootInt, BaseInt, ExpRing> {
    _tp: std::marker::PhantomData<ExpRing>,
    /// The LWE dimension, for binary LWE secrets
    pub(crate) dim_lwe: usize,
    /// log2(N) where the NTRU ciphertext are defined over mod X^N + 1
    pub(crate) log_deg_ntru: u64,
    /// The number of components for gadget decomposition
    pub(crate) dim_ngs: usize,
    /// The base for gadget decomposition, `gadget_base`^`dim_ngs` >= BootInt modulus
    pub(crate) gadget_base: BootInt,
    /// The dimension of the gadget decomposition in the NTRU->LWE KSK
    pub(crate) ksk_ntru_lwe_dim: usize,
    /// The base for the decomposition used in the NTRU->LWE KSK
    pub(crate) ksk_ntru_lwe_base: BootInt,
    /// The (discrete) distribution used to sample LWE noise
    pub(crate) err_stdev_lwe: f64,
    /// The (discrete) distribution used for fresh NTRU noise polynomials
    pub(crate) err_stdev_ntru: f64,
    /// The scale factor to embed the message in an LWE ciphertext
    pub(crate) scale_lwe: BaseInt,
    /// Half the scale factor to embed the message in an LWE ciphertext
    pub(crate) half_scale_lwe: BaseInt,
    /// The scale factor to embed the message in an NGS ciphertext
    pub(crate) scale_ntru: BootInt,
    /// The scale factor for the NTRU key: f = f' * `scale_ntru_key` + 1
    /// And as such, `scale_ntru` * `scale_ntru_key` ≈ BootInt modulus
    pub(crate) scale_ntru_key: BootInt,
    /// Half the scale factor to embed the message in an NGS ciphertext
    pub(crate) half_scale_ntru: BootInt,

    // NOTE: Putting the lazylock on this doesn't entirely solve
    // the issue of caching initialization, as (in theory) you could
    // use multiple `Params` instances with the same ntru degree and
    // hence duplicate FFTPlan work/space.
    // It's fine in practice though :)
    /// The fft plan/engine to use
    pub(crate) fft: crate::fft::FFTPlan,
}

/// Parameters that define the size, security and efficiency of the entire FINAL FHE scheme
pub type Params<BootInt, BaseInt, ExpRing> = LazyLock<ParamsInner<BootInt, BaseInt, ExpRing>>;

#[cfg(test)]
pub(crate) type TESTTYPE = M31;

// TODO: check!
/// Some parameters to test the scheme with during development
#[cfg(test)]
pub static TESTPARAMS: Params<TESTTYPE, TESTTYPE, Z2k<12>> = Params::new(|| ParamsInner {
    _tp: std::marker::PhantomData,
    dim_lwe: 930,
    log_deg_ntru: 10,
    dim_ngs: 3,
    gadget_base: <M31 as crate::modular::Modular>::new_unchecked(1 << 11),
    ksk_ntru_lwe_dim: 3,
    ksk_ntru_lwe_base: <M31 as crate::modular::Modular>::new_unchecked(1 << 11),
    err_stdev_lwe: 4.39,
    err_stdev_ntru: 4.39,
    scale_lwe: <M31 as crate::modular::Modular>::new_unchecked(536870912),
    half_scale_lwe: <M31 as crate::modular::Modular>::new_unchecked(268435456),
    scale_ntru: <M31 as crate::modular::Modular>::new_unchecked(536870912),
    half_scale_ntru: <M31 as crate::modular::Modular>::new_unchecked(268435456),
    scale_ntru_key: crate::modular::int_to_field(4u8.into()),

    fft: crate::fft::FFTPlan::new(1 << 10),
});

// TODO: check, this is entirely random garbage
pub static LARGEPARAMS: Params<M61, M61, Z2k<16>> = Params::new(|| ParamsInner {
    _tp: std::marker::PhantomData,
    dim_lwe: 2048,
    log_deg_ntru: 14,
    dim_ngs: 4,
    gadget_base: <M61 as crate::modular::Modular>::new_unchecked(16),
    ksk_ntru_lwe_dim: 4,
    ksk_ntru_lwe_base: <M61 as crate::modular::Modular>::new_unchecked(16),
    err_stdev_lwe: 12.0,
    err_stdev_ntru: 12.0,
    scale_lwe: <M61 as crate::modular::Modular>::new_unchecked(576460752303423488),
    half_scale_lwe: <M61 as crate::modular::Modular>::new_unchecked(288230376151711744),
    scale_ntru: <M61 as crate::modular::Modular>::new_unchecked(576460752303423488),
    half_scale_ntru: <M61 as crate::modular::Modular>::new_unchecked(288230376151711744),
    scale_ntru_key: crate::modular::int_to_field(4u8.into()),

    fft: crate::fft::FFTPlan::new(1 << 14),
});

//! Parameters defining the entire scheme

use std::ops::IndexMut;

#[cfg(test)]
use crate::modular::M31;

use crate::modular::{Modular, Z2k, Q17, Q20, Q30, Q35};

/// A cryptographically secure RNG trait for convenience
pub trait Rng: rand::Rng + rand::CryptoRng {}
impl<T: rand::Rng + rand::CryptoRng> Rng for T {}

pub trait Params {
    /// The type used for bootstrapping values (NTRU)
    type BootInt: Modular + std::fmt::Debug;
    /// The type used for the base encryption scheme (LWE)
    type BaseInt: Modular + std::fmt::Debug;
    /// The type use to modswitch for homomorphic decryption during bootstrapping
    type ExpRing: Modular + std::fmt::Debug;

    /// The LWE dimension, for binary LWE secrets
    const DIM_LWE: usize;
    /// log2(N) where the NTRU ciphertext are defined over mod X^N + 1
    const LOG_DEG_NTRU: usize;
    /// The number of components for gadget decomposition
    const DIM_NGS: usize;
    /// The base for gadget decomposition, `gadget_base`^`dim_ngs` >= BootInt modulus
    fn gadget_base() -> Self::BootInt;
    /// The dimension of the gadget decomposition in the NTRU->LWE KSK
    const KSK_NTRU_LWE_DIM: usize;
    /// The base for the decomposition used in the NTRU->LWE KSK
    fn ksk_ntru_lwe_base() -> Self::BootInt;
    /// The dimension of the gadget decomposition in the LWE->LWE KSK
    const KSK_LWE_LWE_DIM: usize;
    /// The base for the decomposition used in the LWE->LWE KSK
    fn ksk_lwe_lwe_base() -> Self::BaseInt;
    /// The (discrete) distribution used to sample LWE noise
    const ERR_STDEV_LWE: f64;
    /// The (discrete) distribution used for fresh NTRU noise polynomials
    const ERR_STDEV_NTRU: f64;
    /// The scale factor to embed the message in an LWE ciphertext
    fn scale_lwe() -> Self::BaseInt;
    /// Half the scale factor to embed the message in an LWE ciphertext
    fn half_scale_lwe() -> Self::BaseInt;
    /// The scale factor to embed the message in an NGS ciphertext
    fn scale_ntru() -> Self::BootInt;
    /// Half the scale factor to embed the message in an NGS ciphertext
    fn half_scale_ntru() -> Self::BootInt;
    /// The scale factor for the NTRU key: f = f' * `scale_ntru_key` + 1
    /// And as such, `scale_ntru` * `scale_ntru_key` ≈ BootInt modulus
    fn scale_ntru_key() -> Self::BootInt;

    /// Get the FFT plan for this size NTRU ciphertexts
    /// Acquiring requires some synchronization, so less is more
    #[inline]
    fn fft() -> &'static crate::fft::FFTPlan {
        use crate::fft::FFTPlan;
        use std::cell::UnsafeCell;
        use std::sync::Mutex;

        // SAFETY: we leak the memory, but it should always be an FFT plan of the right size
        // We do some funky things
        debug_assert!(
            Self::LOG_DEG_NTRU < 20,
            "Add more pointer space for FFT plans"
        );
        // Store as usize, because pointers aren't Sync
        static BY_SIZE: Mutex<UnsafeCell<[usize; 20]>> = Mutex::new(UnsafeCell::new([0; 20]));
        let mut array = BY_SIZE.lock().expect("Mutex poisoned");
        let ptr_storage_ref = array.get_mut().index_mut(Self::LOG_DEG_NTRU);
        if *ptr_storage_ref == 0 {
            *ptr_storage_ref = Box::leak(Box::new(FFTPlan::new(1 << Self::LOG_DEG_NTRU)))
                as *const FFTPlan as usize;
        }
        unsafe { &*(*ptr_storage_ref as *const FFTPlan) }
    }
}

#[cfg(test)]
pub enum WeirdParams {}
#[cfg(test)]
impl Params for WeirdParams {
    type BootInt = M31;
    type BaseInt = M31;
    type ExpRing = Z2k<12>;

    const DIM_LWE: usize = 500;
    const LOG_DEG_NTRU: usize = 10;
    const DIM_NGS: usize = 3;
    #[inline]
    fn gadget_base() -> Self::BootInt {
        Self::BootInt::new_unchecked(1 << 11)
    }
    const KSK_NTRU_LWE_DIM: usize = 20;
    #[inline]
    fn ksk_ntru_lwe_base() -> Self::BootInt {
        Self::BootInt::new_unchecked(3)
    }
    const KSK_LWE_LWE_DIM: usize = 3;
    #[inline]
    fn ksk_lwe_lwe_base() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(1 << 11)
    }
    const ERR_STDEV_LWE: f64 = 4.39;
    const ERR_STDEV_NTRU: f64 = 4.39;
    #[inline]
    fn scale_lwe() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(536870912)
    }
    #[inline]
    fn half_scale_lwe() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(268435456)
    }
    #[inline]
    fn scale_ntru() -> Self::BootInt {
        Self::BootInt::new_unchecked(536870912)
    }
    #[inline]
    fn half_scale_ntru() -> Self::BootInt {
        Self::BootInt::new_unchecked(268435456)
    }
    #[inline]
    fn scale_ntru_key() -> Self::BootInt {
        Self::BootInt::new_unchecked(4)
    }
}

pub enum FinalParams {}
impl Params for FinalParams {
    type BootInt = Q20;
    type BaseInt = Q17;
    type ExpRing = Z2k<11>;

    const DIM_LWE: usize = 610;
    const LOG_DEG_NTRU: usize = 10;
    const DIM_NGS: usize = 5;
    #[inline]
    fn gadget_base() -> Self::BootInt {
        Self::BootInt::new_unchecked(16)
    }
    const KSK_NTRU_LWE_DIM: usize = 13;
    #[inline]
    fn ksk_ntru_lwe_base() -> Self::BootInt {
        Self::BootInt::new_unchecked(3)
    }
    const KSK_LWE_LWE_DIM: usize = 6;
    #[inline]
    fn ksk_lwe_lwe_base() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(8)
    }
    const ERR_STDEV_LWE: f64 = 4.39;
    const ERR_STDEV_NTRU: f64 = 4.39;
    #[inline]
    fn scale_lwe() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(23171)
    }
    #[inline]
    fn half_scale_lwe() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(11585)
    }
    #[inline]
    fn scale_ntru() -> Self::BootInt {
        Self::BootInt::new_unchecked(228207)
    }
    #[inline]
    fn half_scale_ntru() -> Self::BootInt {
        Self::BootInt::new_unchecked(114104)
    }
    #[inline]
    fn scale_ntru_key() -> Self::BootInt {
        Self::BootInt::new_unchecked(4)
    }
}

pub enum MKFinalParams {}
impl Params for MKFinalParams {
    type BootInt = Q30;
    type BaseInt = Q35;
    type ExpRing = Z2k<12>;

    const DIM_LWE: usize = 660;
    const LOG_DEG_NTRU: usize = 11;
    const DIM_NGS: usize = 8;
    #[inline]
    fn gadget_base() -> Self::BootInt {
        Self::BootInt::new_unchecked(16)
    }
    const KSK_NTRU_LWE_DIM: usize = 23;
    #[inline]
    fn ksk_ntru_lwe_base() -> Self::BootInt {
        Self::BootInt::new_unchecked(3)
    }
    const KSK_LWE_LWE_DIM: usize = 9;
    #[inline]
    fn ksk_lwe_lwe_base() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(16)
    }
    const ERR_STDEV_LWE: f64 = 4.0;
    const ERR_STDEV_NTRU: f64 = 0.0;
    #[inline]
    fn scale_lwe() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(8589934605)
    }
    #[inline]
    fn half_scale_lwe() -> Self::BaseInt {
        Self::BaseInt::new_unchecked(4294967303)
    }
    #[inline]
    fn scale_ntru() -> Self::BootInt {
        Self::BootInt::new_unchecked(268435457)
    }
    #[inline]
    fn half_scale_ntru() -> Self::BootInt {
        Self::BootInt::new_unchecked(134217728)
    }
    #[inline]
    fn scale_ntru_key() -> Self::BootInt {
        Self::BootInt::new_unchecked(4)
    }
}

#[cfg(test)]
pub use MKFinalParams as TestParams;
#[cfg(test)]
pub type TESTTYPE = M31;

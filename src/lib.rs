// TODO: this might help moving away from the single-limb assumption
// should it ever be needed, since then you can do a
// `x.into_int::<{ F::MIN_LIMBS_NEEDED }>()`
//
// #![feature(generic_const_exprs)]
//
// Alternatively, to avoid nightly rust,
// we'd need to make a const LIMBS: usize part of the generics of Params

pub mod ciphertext;
pub mod key;
pub mod multiparty;
pub mod params;

pub(crate) mod fft;
pub(crate) mod modular;
pub(crate) mod poly;
#[cfg(test)]
pub(crate) mod test_utils;

#[cfg(feature = "bench")]
pub mod bench {
    pub mod fft {
        pub use crate::fft::*;
    }

    pub mod modular {
        pub use crate::modular::*;
    }

    pub mod poly {
        pub use crate::poly::*;
    }
}

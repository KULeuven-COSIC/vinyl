// Note, if things turn out to be too slow, maybe switch everything to flint nmod + nmod_poly
use std::ops::{Div, Sub};

use num_modular::{ModularInteger, ModularUnaryOps};

pub trait ModSwitch<Target> {
    fn modswitch(self) -> Target;
}

impl<F: Modular, T: Modular> ModSwitch<T> for F {
    fn modswitch(self) -> T {
        todo!()
    }
}

pub trait UnaryOps {
    fn inv(&self) -> Self;
}

pub trait Modular: num_modular::ModularInteger + UnaryOps + Clone {
    fn construct(n: Self::Base) -> Self;
    fn modulus() -> Self::Base;

    fn centered(&self) -> Self::Base
    where
        Self::Base: Div<usize, Output = Self::Base>
            + Sub<Self::Base, Output = Self::Base>
            + PartialOrd
            + Copy,
    {
        let m = self.modulus();
        let r = self.residue();
        if r > m / 2 {
            r - m
        } else {
            r
        }
    }
}

#[derive(Debug)]
pub struct Params<BootInt, BaseInt, NTRUExpInt> {
    _tp: std::marker::PhantomData<(BootInt, BaseInt, NTRUExpInt)>,
    pub(crate) dim_lwe: u64,
    pub(crate) log_deg_ntru: u64,
    pub(crate) dim_ngs: u64,

    pub(crate) mod_plain: u64,
    pub(crate) mod_carry: u64,

    pub(crate) err_stdev: f64,
}

pub type Q64 = num_modular::FixedMersenneInt<64, 59>;
type Mod2PowN<const N: u8> = num_modular::FixedMersenneInt<N, 0>;

impl<const P: u8, const K: u128> UnaryOps for num_modular::FixedMersenneInt<P, K> {
    fn inv(&self) -> Self {
        let m = self.modulus();
        let divisor = self.residue().invm(&m).expect("Modular inverse failure");
        self.convert(divisor)
    }
}

impl<const P: u8, const K: u128> Modular for num_modular::FixedMersenneInt<P, K> {
    fn construct(n: Self::Base) -> Self {
        Self::new(n, &num_modular::FixedMersenne::<P, K>::MODULUS)
    }

    fn modulus() -> Self::Base {
        num_modular::FixedMersenne::<P, K>::MODULUS
    }
}

// TODO: check security and optimality :)
pub const MESSAGE_7_CARRY_0: Params<Q64, Q64, Mod2PowN<15>> = Params {
    _tp: std::marker::PhantomData,
    dim_lwe: 930,
    log_deg_ntru: 14,
    dim_ngs: 3,

    mod_plain: 128,
    mod_carry: 2,

    err_stdev: 4.39,
};

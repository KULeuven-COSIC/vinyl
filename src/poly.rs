use std::ops::{Add, Div, Mul, Sub};

use crate::modular::Modular;
use concrete_fft::fft128::f128;

// NOTE: Might want to switch to FLINT if too inefficient?
/// A polynomial, stored as a coefficient vector, constant term first
#[derive(Debug, Clone)]
pub struct Poly<T>(Vec<T>);

impl<T> Poly<T>
where
    T: std::fmt::Debug,
    T: Modular,
    T::Base: From<u64> + TryInto<u64>,
    T::Base: PartialOrd,
{
    /// Return the `inv: Poly<T>` such that `inv * self = 1 (mod x^n + 1)`
    /// None when not invertible
    pub fn invert(&self) -> Option<Poly<T>> {
        use flint_sys::deps::mp_limb_t;
        use flint_sys::nmod_poly::*;

        debug_assert!(self.0.len() > 0);
        debug_assert!(self.0[0].modulus() <= From::from(mp_limb_t::MAX));
        debug_assert!(mp_limb_t::BITS >= u64::BITS);

        let modulus: mp_limb_t = self.0[0]
            .modulus()
            .try_into()
            .map_err(|_| ())
            .expect("Should be caught by debug_assert");
        let mut me = unsafe {
            let mut me = std::mem::MaybeUninit::uninit();
            nmod_poly_init2(me.as_mut_ptr(), modulus, self.0.len() as i64);
            for (i, x) in self.iter().enumerate() {
                nmod_poly_set_coeff_ui(
                    me.as_mut_ptr(),
                    i as i64,
                    x.residue().try_into().map_err(|_| ()).unwrap(),
                );
            }
            me.assume_init()
        };
        let mut defpoly = unsafe {
            let mut r = std::mem::MaybeUninit::uninit();
            nmod_poly_init2(r.as_mut_ptr(), modulus, self.0.len() as i64 + 1);
            nmod_poly_set_coeff_ui(r.as_mut_ptr(), 0, 1);
            nmod_poly_set_coeff_ui(r.as_mut_ptr(), self.0.len() as i64, 1);
            r.assume_init()
        };
        let mut res = unsafe {
            let mut r = std::mem::MaybeUninit::uninit();
            nmod_poly_init(r.as_mut_ptr(), modulus);
            r.assume_init()
        };
        let status = unsafe {
            nmod_poly_invmod(
                &mut res as *mut _,
                &mut me as *mut _,
                &mut defpoly as *mut _,
            )
        };

        let inv = if status == 0 {
            None
        } else {
            let mut v = Vec::with_capacity(self.0.len());
            for i in 0..self.0.len() {
                v.push(T::construct(From::from(unsafe {
                    nmod_poly_get_coeff_ui(&mut res as *mut _, i as i64)
                })));
            }
            Some(Self(v))
        };

        unsafe {
            nmod_poly_clear(&mut me as *mut _);
            nmod_poly_clear(&mut defpoly as *mut _);
            nmod_poly_clear(&mut res as *mut _);
        }
        inv
    }
}

impl<T> Poly<T>
where
    T: Modular,
    T::Base: From<u64>,
{
    /// Create a new polynomial of degree `n`, initialized to the all-zero polynomial
    pub fn new(n: u64) -> Self {
        Self(vec![T::construct(From::from(0)); n as usize])
    }

    /// Iterate over coefficients, constant term first
    pub fn iter(&self) -> impl Iterator<Item = &T> + '_ {
        self.0.iter()
    }

    /// Iterate over coefficients, constant term first; mutably
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> + '_ {
        self.0.iter_mut()
    }
}

impl<T> IntoIterator for Poly<T> {
    type Item = T;

    type IntoIter = <Vec<T> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

/// Allow to define traits that act pointwise (or scalar) over the elements of something iterable
macro_rules! pointwise(
    ($trt: ident, $fn: ident, $t: ident) => {
        impl<T> $trt<$t<T>> for $t<T>
        where
            T: $trt<T, Output=T>,
        {
            type Output = Self;
            fn $fn(self, other: Self) -> Self {
                debug_assert!(self.0.len() == other.0.len());
                Self(
                    self.into_iter()
                        .zip(other.into_iter())
                        .map(|(a, b)| $trt::$fn(a, b))
                        .collect()
                )
            }
        }
    };
    (scalar, $trt: ident, $fn: ident, $t: ident) => {
        impl<T> $trt<T> for $t<T>
        where
            T: $trt<T, Output=T> + Clone,
        {
            type Output = Self;
            fn $fn(self, other: T) -> Self {
                Self(
                    self.into_iter()
                        .map(|a| $trt::$fn(a, other.clone()))
                        .collect()
                )
            }
        }
    }
);

pointwise!(Add, add, Poly);
pointwise!(Sub, sub, Poly);
pointwise!(scalar, Mul, mul, Poly);
pointwise!(scalar, Div, div, Poly);

/// A polynomial stored in FFT form
#[derive(Debug, Clone)]
pub struct FFTPoly(Vec<f128>);

use std::ops::{Add, Div, Mul, Sub};

use swanky_field::PrimeFiniteField;
use tfhe_fft::fft128::f128;

use crate::params::Rng;

// TODO? Move somewhere dedicated
/// Complex numbers, generic over the base type
/// We implement this ourselves because tfhe_fft::fft128::f128 doesn't implement `Num`,
/// but `num_complex::Complex<T>` requires it to provide useful operations
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Complex<T> {
    pub re: T,
    pub im: T,
}

impl<T> Complex<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }
}

impl<T> Add<&Complex<T>> for Complex<T>
where
    T: for<'a> Add<&'a T, Output = T>,
{
    type Output = Complex<T>;

    fn add(self, rhs: &Complex<T>) -> Self::Output {
        Self {
            re: self.re + &rhs.re,
            im: self.im + &rhs.im,
        }
    }
}

impl<T> Add<Complex<T>> for Complex<T>
where
    T: Add<T, Output = T>,
{
    type Output = Complex<T>;

    fn add(self, rhs: Complex<T>) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

// Kinda sucks that we need to clone more, and can't mult f128 * &f128; but oh well
impl<T> Mul<&Complex<T>> for Complex<T>
where
    T: Mul<T, Output = T> + Sub<T, Output = T> + Add<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    fn mul(self, rhs: &Complex<T>) -> Self::Output {
        Self {
            re: self.re.clone() * rhs.re.clone() - self.im.clone() * rhs.im.clone(),
            im: self.re * rhs.im.clone() + self.im * rhs.re.clone(),
        }
    }
}

impl<T> Mul<Complex<T>> for Complex<T>
where
    T: Mul<T, Output = T> + Sub<T, Output = T> + Add<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    fn mul(self, rhs: Complex<T>) -> Self::Output {
        Self {
            re: self.re.clone() * rhs.re.clone() - self.im.clone() * rhs.im.clone(),
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl<T> Mul<&T> for Complex<T>
where
    T: for<'a> Mul<&'a T, Output = T>,
{
    type Output = Complex<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl<T> Mul<T> for Complex<T>
where
    T: Mul<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self {
            re: self.re * rhs.clone(),
            im: self.im * rhs,
        }
    }
}

// NOTE: Might want to switch to FLINT if too inefficient?
/// A polynomial, stored as a coefficient vector, constant term first
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly<T>(pub(crate) Vec<T>);

impl<F> Poly<F>
where
    F: PrimeFiniteField,
{
    const _LIMB_NUM: () = assert!(F::MIN_LIMBS_NEEDED == 1);
    const _LIMB_SIZE: () =
        assert!(crypto_bigint::Limb::BITS <= flint_sys::deps::mp_limb_t::BITS as usize);

    /// Return the `inv: Poly<T>` such that `inv * self = 1 (mod x^n + 1)`
    /// None when not invertible
    pub fn invert(&self) -> Option<Poly<F>> {
        use flint_sys::deps::mp_limb_t;
        use flint_sys::nmod_poly::*;

        debug_assert!(!self.0.is_empty());

        let modulus: mp_limb_t = F::modulus_int().into();

        let mut me = unsafe {
            let mut me = std::mem::MaybeUninit::uninit();
            nmod_poly_init2(me.as_mut_ptr(), modulus, self.0.len() as i64);
            for (i, x) in self.iter().enumerate() {
                nmod_poly_set_coeff_ui(me.as_mut_ptr(), i as i64, x.into_int().into());
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
                let coef = unsafe { nmod_poly_get_coeff_ui(&mut res as *mut _, i as i64) };
                v.push(F::try_from_int::<1>(coef.into()).expect("Non-reduced element from flint"));
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

impl<T> Poly<T> {
    /// Iterate over coefficients, constant term first
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.0.iter()
    }

    /// Iterate over coefficients, constant term first; mutably
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> + '_ {
        self.0.iter_mut()
    }

    /// What's the degree of this polynomial?
    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }
}

impl<T> Poly<T>
where
    T: swanky_field::FiniteRing,
{
    /// Create a new polynomial of degree `n - 1`, initialized to the all-zero polynomial
    pub fn new(n: usize) -> Self {
        Self(vec![T::ZERO; n])
    }

    /// Generate a random polynomial over `T` of degree `n - 1`
    pub fn rand(n: usize, rng: &mut impl Rng) -> Self {
        Self(Vec::from_iter(
            std::iter::from_fn(|| Some(T::random(rng))).take(n),
        ))
    }

    /// Get a specific coefficient
    pub(crate) fn get(&self, i: usize) -> T {
        self.0[i]
    }
}

impl<T> IntoIterator for Poly<T> {
    type Item = T;

    type IntoIter = <Vec<T> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Poly<T> {
    type Item = &'a T;

    type IntoIter = <&'a [T] as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        (&self.0).into_iter()
    }
}

impl<T> std::iter::Sum for Poly<T>
where
    T: Add<T, Output = T> + Default + Clone,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut pkbl = iter.peekable();
        let Some(first) = pkbl.peek() else {
            return Self(vec![]);
        };
        let mut res = Self(vec![T::default(); first.0.len()]);

        for el in pkbl {
            res = res + el;
        }

        res
    }
}

/// Allow to define traits that act pointwise (or scalar) over the elements of something iterable
macro_rules! pointwise(
    ($trt: ident, $fn: ident, $t: ident) => {
        pointwise!(@impl, $trt, $fn, $t, $t<T>, $trt<T, Output=T>);
        pointwise!(@impl, $trt, $fn, $t, &$t<T>, for<'a> $trt<&'a T, Output=T>);
    };
    (scalar, $trt: ident, $fn: ident, $t: ident) => {
        pointwise!(@impl scalar, $trt, $fn, $t, T, $trt<T, Output=T> + Copy);
        pointwise!(@impl scalar, $trt, $fn, $t, &T, for<'a> $trt<&'a T, Output=T>);
    };
    (@impl, $trt: ident, $fn: ident, $t: ident, $target:ty, $($bound:tt)*) => {
        impl<T> $trt<$target> for $t<T>
        where
            T: $($bound)*,
        {
            type Output = Self;
            fn $fn(self, other: $target) -> Self {
                debug_assert_eq!(self.0.len(), other.0.len());
                Self(
                    self.into_iter()
                        .zip(other.into_iter())
                        .map(|(a, b)| $trt::$fn(a, b))
                        .collect()
                )
            }
        }
    };
    (@impl scalar, $trt: ident, $fn: ident, $t: ident, $target:ty, $($bound:tt)*) => {
        impl<T> $trt<$target> for $t<T>
        where
            T: $($bound)*
        {
            type Output = Self;
            fn $fn(self, other: $target) -> Self {
                Self(
                    self.into_iter()
                        .map(|a| $trt::$fn(a, other))
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

impl<F: PrimeFiniteField> Add<F> for Poly<F> {
    type Output = Poly<F>;

    fn add(mut self, rhs: F) -> Self::Output {
        self.0[0] += rhs;
        self
    }
}

/// A polynomial stored in FFT form
#[derive(Debug, Clone)]
pub struct FFTPolyGeneric<T>(pub(crate) Vec<T>);
pub type FFTPoly = FFTPolyGeneric<Complex<f128>>;

impl<T> IntoIterator for FFTPolyGeneric<T> {
    type Item = T;

    type IntoIter = <Vec<T> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a FFTPolyGeneric<T> {
    type Item = &'a T;

    type IntoIter = <&'a [T] as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        (&self.0).into_iter()
    }
}

impl<T> FFTPolyGeneric<T> {
    pub fn new(x: Vec<T>) -> Self {
        Self(x)
    }
}

impl std::iter::Sum for FFTPoly {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut pkbl = iter.peekable();
        let Some(first) = pkbl.peek() else {
            return Self(vec![]);
        };
        let mut res = Self(vec![
            Complex {
                re: f128(0.0, 0.0),
                im: f128(0.0, 0.0)
            };
            first.0.len()
        ]);

        for el in pkbl {
            res = res + el;
        }

        res
    }
}

pointwise!(Mul, mul, FFTPolyGeneric);
pointwise!(Add, add, FFTPolyGeneric);
pointwise!(scalar, Mul, mul, FFTPolyGeneric);
pointwise!(scalar, Add, add, FFTPolyGeneric);

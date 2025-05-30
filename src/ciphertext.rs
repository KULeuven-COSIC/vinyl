#[cfg(test)]
use swanky_field::FiniteRing;
use swanky_field::PrimeFiniteField;

use crate::{
    modular::{int_to_field, ModSwitch},
    params::Params,
    poly::{FFTPoly, Poly},
};

#[derive(Clone, Debug)]
pub struct MKLweCiphertext<F, const N: usize> {
    pub(crate) a: [Vec<F>; N],
    pub(crate) b: F,
}

// TODO: implement some "claim" functionality or smth to get the single sample out of the array
pub type LweCiphertext<F> = MKLweCiphertext<F, 1>;

impl<F> LweCiphertext<F> {
    /// Obtain the single `(a, b)` sample from this ciphertext.
    /// Useful to avoid pain with the `a` being stored in an array
    pub(crate) fn unpack(self) -> (Vec<F>, F) {
        // SAFETY: Single element array, layout works
        (
            unsafe { std::mem::transmute::<[Vec<F>; 1], Vec<F>>(self.a) },
            self.b,
        )
    }
}

impl<F, const N: usize> MKLweCiphertext<F, N>
where
    F: PrimeFiniteField,
{
    /// Decrypt a ciphertext with explicit parameters, rather than reading from a `Params`
    /// Useful e.g. when dealing with a different field
    pub(crate) fn decrypt_explicit<Key: std::borrow::Borrow<crate::key::LWEKey>>(
        self,
        keys: &[Key; N],
        dim: usize,
        half_scale: F,
        scale: F,
    ) -> u8 {
        #[cfg(debug_assertions)]
        {
            for a in self.a.iter() {
                debug_assert_eq!(a.len(), dim);
            }
            for key in keys {
                debug_assert_eq!(dim, key.borrow().dim);
                debug_assert!(key.borrow().key.len() * crate::key::LweKeyEl::BITS as usize >= dim);
            }
        }
        let mask: F = self
            .a
            .into_iter()
            .zip(keys.iter())
            .flat_map(|(a, key)| {
                a.into_iter()
                    .zip(key.borrow().iter(dim))
                    .map(|(ai, si): (F, F)| ai * si)
            })
            .sum();

        // Add scale / 2 and floor div to round
        let out = (self.b + mask + half_scale).into_int::<1>()
            / crypto_bigint::NonZero::new(scale.into_int::<1>()).unwrap();
        out.bit_vartime(0) as u8
    }

    pub(crate) fn double(self) -> Self {
        let b = self.b + self.b;
        MKLweCiphertext {
            a: self.a.map(|aa| aa.into_iter().map(|ai| ai + ai).collect()),
            b,
        }
    }

    #[allow(private_bounds)]
    pub fn decrypt<P: Params<BaseInt = F>, Key: std::borrow::Borrow<crate::key::LWEKey>>(
        self,
        keys: &[Key; N],
    ) -> u8
    where
        P::BaseInt: PrimeFiniteField,
    {
        self.decrypt_explicit(keys, P::DIM_LWE, P::half_scale_lwe(), P::scale_lwe())
    }
}

impl<F, T, const N: usize> ModSwitch<MKLweCiphertext<T, N>> for MKLweCiphertext<F, N>
where
    F: ModSwitch<T>,
{
    fn modswitch(self) -> MKLweCiphertext<T, N> {
        MKLweCiphertext {
            // Via Poly because we can't impl it for Vec easily
            a: self.a.map(|a| Poly(a).modswitch().0),
            b: self.b.modswitch(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct NtruScalarCiphertext<T> {
    pub(crate) ct: Poly<T>,
}

// TODO: approximate decomposition
#[derive(Clone, Debug)]
pub struct NtruVectorCiphertext {
    pub(crate) ct: Vec<FFTPoly>,
}

impl<T: swanky_field::PrimeFiniteField> NtruScalarCiphertext<T> {
    // TODO: again assumes everything fits into a single limb, also not constant time
    /// Gadget decomposition into `n` parts. Returns an a vec such that
    /// `sum(result[i] * B^i) = self`
    /// If base is known to be a power of two, use `gadget_decomp_pow2` instead
    pub(crate) fn gadget_decomp_generic(self, dim: usize, base: T) -> Vec<Poly<T>> {
        let mut res = Vec::with_capacity(dim);
        let mut as_ints: Vec<_> = self
            .ct
            .0
            .into_iter()
            .map(|x| x.into_int::<1>().as_words()[0])
            .collect();

        let base_int = base.into_int::<1>().as_words()[0];
        debug_assert!(base_int.pow(dim as u32) > T::modulus_int::<1>().as_words()[0]);

        for _ in 0..dim {
            let mut part = Vec::<T>::with_capacity(as_ints.len());
            for coeff in as_ints.iter_mut() {
                part.push(int_to_field((*coeff % base_int).into()));
                *coeff /= base_int;
            }
            res.push(Poly(part));
        }
        debug_assert!(as_ints.into_iter().all(|x| x == 0));

        res
    }

    // TODO: again assumes everything fits into a single limb
    /// Gadget decomposition into `n` parts. Returns an a vec such that
    /// `sum(result[i] * B^i) = self`
    /// Assumes `base` is a power of two!
    pub(crate) fn gadget_decomp_pow2(self, dim: usize, base: T) -> Vec<Poly<T>> {
        let mut res = Vec::with_capacity(dim);
        let mut as_ints: Vec<_> = self
            .ct
            .0
            .into_iter()
            .map(|x| x.into_int::<1>().as_words()[0])
            .collect();

        let base_int = base.into_int::<1>().as_words()[0];
        debug_assert_eq!(base_int.count_ones(), 1);
        debug_assert!(base_int.pow(dim as u32) > T::modulus_int::<1>().as_words()[0]);
        let base_mask = base_int - 1;
        let base_shift = base_int.ilog2();

        for _ in 0..dim {
            let mut part = Vec::<T>::with_capacity(as_ints.len());

            for coeff in as_ints.iter_mut() {
                part.push(
                    T::try_from_int::<1>((*coeff & base_mask).into())
                        .expect("Gadget decomp base guaranteed to fit in modulus"),
                );
                *coeff >>= base_shift;
            }

            res.push(Poly(part));
        }
        res
    }

    /// Compute the external product with a vector encryption
    /// Not implemented as a Mul trait for two reasons:
    /// - We need access to the params
    /// - This is costly, so explicit is better than implicit here
    #[cfg_vis::cfg_vis(feature = "bench", pub)]
    pub(crate) fn external_product<P: Params<BootInt = T>>(
        self,
        rhs: &NtruVectorCiphertext,
        fft: &crate::fft::FFTPlan,
    ) -> Self {
        Self {
            ct: fft.inv(
                self.gadget_decomp_pow2(P::DIM_NGS, P::gadget_base())
                    .into_iter()
                    .zip(rhs.ct.iter())
                    .map(|(x, y)| fft.fwd(x) * y)
                    .sum(),
            ),
        }
    }

    /// Rotate self by `n` places, negacyclically
    pub(crate) fn rot(self, n: usize) -> Self {
        let len = self.ct.0.len();
        let mut n = n % (2 * len);
        // TODO: not constant time
        let sign = if n >= len {
            n -= len;
            -T::ONE
        } else {
            T::ONE
        };
        let mut res = Poly::new(self.ct.0.len());

        for i in 0..(len - n) {
            res.0[i + n] = sign * self.ct.0[i];
        }
        for i in (len - n)..len {
            res.0[i + n - len] = -sign * self.ct.0[i];
        }

        Self { ct: res }
    }
}

impl<T: swanky_field::PrimeFiniteField + for<'a> std::ops::Mul<&'a T, Output = T>>
    NtruScalarCiphertext<T>
{
    #[cfg(any(feature = "bench", test))]
    #[cfg_vis::cfg_vis(feature = "bench", pub)]
    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the "normal" scale $\Delta = Q/4$
    pub(crate) fn trivial<P: Params<BootInt = T>>(pt: Poly<T>) -> Self {
        Self::trivial_scale(pt, &P::scale_ntru())
    }

    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the "halved" scale $\Delta = Q/8$ as used for the test vector
    pub(crate) fn trivial_half<P: Params<BootInt = T>>(pt: Poly<T>) -> Self {
        Self::trivial_scale(pt, &P::half_scale_ntru())
    }

    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the provided scale $\Delta$
    fn trivial_scale(pt: Poly<T>, scale: &T) -> Self {
        Self { ct: pt * scale }
    }
}

impl<T: PrimeFiniteField> std::ops::Add<NtruScalarCiphertext<T>> for NtruScalarCiphertext<T> {
    type Output = Self;

    fn add(self, rhs: NtruScalarCiphertext<T>) -> Self::Output {
        Self {
            ct: self.ct + rhs.ct,
        }
    }
}

impl<T: PrimeFiniteField> std::ops::Sub<NtruScalarCiphertext<T>> for NtruScalarCiphertext<T> {
    type Output = Self;

    fn sub(self, rhs: NtruScalarCiphertext<T>) -> Self::Output {
        Self {
            ct: self.ct - rhs.ct,
        }
    }
}

impl NtruVectorCiphertext {
    #[cfg(test)]
    pub(crate) fn trivial<P: Params>(m: Poly<P::BootInt>) -> Self
    where
        P::BootInt: PrimeFiniteField,
    {
        let mut ct = Vec::with_capacity(P::DIM_NGS);

        let mut base = P::BootInt::ONE;
        for _ in 0..P::DIM_NGS {
            ct.push(P::fft().fwd(m.clone() * base));
            base *= P::gadget_base();
        }

        Self { ct }
    }

    // pub(crate) fn monomial<P: Params>(exponent: usize) -> Self
    // where
    //     P::BootInt: PrimeFiniteField,
    // {
    //     #[allow(non_snake_case)]
    //     let N = 1 << P::LOG_DEG_NTRU;
    //     let mut poly = Poly::new(N);
    //     let exponent = exponent % (2 * N);
    //     let (exponent, value) = if exponent >= N {
    //         (exponent - N, -P::BootInt::ONE)
    //     } else {
    //         (exponent, P::BootInt::ONE)
    //     };
    //     poly.0[exponent] = value;
    //     Self::trivial::<P>(poly)
    // }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::params::{Params, TestParams};
    use crate::test_utils::*;

    type BootInt = <TestParams as Params>::BootInt;

    use swanky_field::FiniteRing;

    #[test]
    fn gadget_decomp_consistency_generic() {
        let rng = &mut rng(None);
        let poly = Poly::<BootInt>::rand(1 << TestParams::LOG_DEG_NTRU, rng);
        let decomp = NtruScalarCiphertext { ct: poly.clone() }.gadget_decomp_generic(
            TestParams::KSK_NTRU_LWE_DIM,
            TestParams::ksk_ntru_lwe_base(),
        );
        assert_eq!(decomp.len(), TestParams::KSK_NTRU_LWE_DIM);
        let sum = decomp
            .into_iter()
            .enumerate()
            .map(|(i, x)| x * TestParams::ksk_ntru_lwe_base().pow_var_time(i as u128))
            .sum::<Poly<_>>();
        assert_eq!(poly, sum);
    }

    #[test]
    fn gadget_decomp_consistency_pow2() {
        let rng = &mut rng(None);
        let poly = Poly::<BootInt>::rand(1 << TestParams::LOG_DEG_NTRU, rng);
        let decomp = NtruScalarCiphertext { ct: poly.clone() }
            .gadget_decomp_pow2(TestParams::DIM_NGS, TestParams::gadget_base());
        assert_eq!(decomp.len(), TestParams::DIM_NGS);
        let sum = decomp
            .into_iter()
            .enumerate()
            .map(|(i, x)| x * TestParams::gadget_base().pow_var_time(i as u128))
            .sum::<Poly<_>>();
        assert_eq!(poly, sum);
    }

    #[test]
    fn external_product_trivial() {
        let rng = &mut rng(None);

        for _ in 0..5 {
            let scalar = Poly::<BootInt>::rand(1 << TestParams::LOG_DEG_NTRU, rng);
            let x = NtruScalarCiphertext::trivial::<TestParams>(scalar.clone());

            let vector = Poly::<BootInt>::rand(1 << TestParams::LOG_DEG_NTRU, rng);
            let y = NtruVectorCiphertext::trivial::<TestParams>(vector.clone());

            let z = x.external_product::<TestParams>(&y, TestParams::fft());
            let poly_prod = NtruScalarCiphertext::trivial::<TestParams>(Poly(
                slow_negacyclic_mult(&scalar.0, &vector.0, BootInt::ZERO),
            ));
            assert_eq!(z.ct, poly_prod.ct);
        }
    }

    #[cfg(not(feature = "bench"))]
    #[test]
    fn external_product_noisy() {
        use crate::modular::{int_to_field, sample_ternary};
        use rand::Rng;

        let rng = &mut rng(None);
        let len = 1 << TestParams::LOG_DEG_NTRU;
        let (key, _) = crate::key::NTRUKey::new::<TestParams>(rng);
        let fft = TestParams::fft();

        for _ in 0..20 {
            for bit in 0..=1 {
                // accumulator is a ternary plaintext
                let mut scalar = Poly::<BootInt>::new(len);
                scalar.iter_mut().for_each(|x| *x = sample_ternary(rng));
                let x = NtruScalarCiphertext::trivial::<TestParams>(scalar.clone());

                // vector ciphertexts are ternary monomials
                let mut vector1 = Poly::<BootInt>::new(len);
                vector1.0[rng.gen_range(0..len)] = sample_ternary(rng);
                let y = key.enc_vec::<TestParams>(&vector1, rng);

                let mut vector2 = Poly::<BootInt>::new(len);
                vector2.0[rng.gen_range(0..len)] = sample_ternary(rng);
                let z = key.enc_vec::<TestParams>(&vector2, rng);

                let mut bit_vec = vec![BootInt::ZERO; len];
                bit_vec[0] = int_to_field(bit.into());
                let b = key.enc_bit_vec::<TestParams>(bit, rng);

                let res = x
                    .external_product::<TestParams>(&y, fft)
                    .external_product::<TestParams>(&z, fft)
                    .external_product::<TestParams>(&b, fft);
                let dec = key.dec_scalar::<TestParams>(res);
                let poly_prod = slow_negacyclic_mult(
                    &slow_negacyclic_mult(
                        &slow_negacyclic_mult(&scalar.0, &vector1.0, BootInt::ZERO),
                        &vector2.0,
                        BootInt::ZERO,
                    ),
                    &bit_vec,
                    BootInt::ZERO,
                );
                assert_eq!(dec, Poly(poly_prod));
            }
        }
    }
}

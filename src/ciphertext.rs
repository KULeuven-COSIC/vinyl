use swanky_field::PrimeFiniteField;

use crate::{
    modular::ModSwitch,
    poly::{FFTPoly, Poly},
};

#[derive(Clone, Debug)]
pub struct LweCiphertext<F> {
    pub(crate) a: Vec<F>,
    pub(crate) b: F,
}

impl<F, T> ModSwitch<LweCiphertext<T>> for LweCiphertext<F>
where
    F: ModSwitch<T>,
{
    fn modswitch(self) -> LweCiphertext<T> {
        LweCiphertext {
            // Via Poly because we can't impl it for Vec easily
            a: Poly(self.a).modswitch().0,
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
    // TODO: again assumes everything fits into a single limb
    /// Gadget decomposition into `n` parts. Returns an a vec such that
    /// `sum(result[i] * B^i) = self`
    pub(crate) fn gadget_decomp(self, dim: usize, base: T) -> Vec<Poly<T>> {
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

            res.push(Poly(part))
        }
        res
    }

    /// Compute the external product with a vector encryption
    /// Not implemented as a Mul trait for two reasons:
    /// - We need access to the params
    /// - This is costly, so explicit is better than implicit here
    pub(crate) fn external_product<BaI, ER>(
        self,
        rhs: &NtruVectorCiphertext,
        params: &crate::params::Params<T, BaI, ER>,
    ) -> Self {
        Self {
            ct: params.fft.inv(
                self.gadget_decomp(params.dim_ngs, params.gadget_base)
                    .into_iter()
                    .zip(rhs.ct.iter())
                    .map(|(x, y)| params.fft.fwd(x) * y)
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

        for i in 0..len {
            res.0[(i + n) % len] = sign * self.ct.0[i];
        }

        Self { ct: res }
    }
}

impl<T: swanky_field::PrimeFiniteField + for<'a> std::ops::Mul<&'a T, Output = T>>
    NtruScalarCiphertext<T>
{
    #[cfg(test)]
    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the "normal" scale $\Delta = Q/4$
    pub(crate) fn trivial<BaI, ER>(
        pt: Poly<T>,
        params: &crate::params::Params<T, BaI, ER>,
    ) -> Self {
        Self::trivial_scale(pt, &params.scale_ntru)
    }

    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the "halved" scale $\Delta = Q/8$ as used for the test vector
    pub(crate) fn trivial_half<BaI, ER>(
        pt: Poly<T>,
        params: &crate::params::Params<T, BaI, ER>,
    ) -> Self {
        Self::trivial_scale(pt, &params.half_scale_ntru)
    }

    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the provided scale $\Delta$
    fn trivial_scale(pt: Poly<T>, scale: &T) -> Self {
        Self { ct: pt * scale }
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
    pub(crate) fn trivial<T: PrimeFiniteField, BaI, ER>(
        m: Poly<T>,
        params: &crate::params::Params<T, BaI, ER>,
    ) -> Self {
        let mut ct = Vec::with_capacity(params.dim_ngs);

        let mut base = T::ONE;
        for _ in 0..params.dim_ngs {
            ct.push(params.fft.fwd(m.clone() * base));
            base = base * params.gadget_base;
        }

        Self { ct }
    }

    pub(crate) fn monomial<T: PrimeFiniteField, BaI, ER>(
        exponent: usize,
        params: &crate::params::Params<T, BaI, ER>,
    ) -> Self {
        #[allow(non_snake_case)]
        let N = 1 << params.log_deg_ntru;
        let mut poly = Poly::new(N);
        let exponent = exponent % (2 * N);
        let (exponent, value) = if exponent >= N {
            (exponent - N, -T::ONE)
        } else {
            (exponent, T::ONE)
        };
        poly.0[exponent] = value;
        Self::trivial(poly, params)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::modular::{int_to_field, sample_ternary};
    use crate::params::{TESTPARAMS, TESTTYPE};
    use crate::test_utils::*;

    use rand::Rng;
    use swanky_field::FiniteRing;

    #[test]
    fn gadget_decomp_consistency() {
        let rng = &mut rng(None);
        let poly = Poly::<TESTTYPE>::rand(1 << TESTPARAMS.log_deg_ntru, rng);
        let decomp = NtruScalarCiphertext { ct: poly.clone() }
            .gadget_decomp(TESTPARAMS.dim_ngs, TESTPARAMS.gadget_base);
        let sum = decomp
            .into_iter()
            .enumerate()
            .map(|(i, x)| x * TESTPARAMS.gadget_base.pow_var_time(i as u128))
            .sum::<Poly<_>>();
        assert_eq!(poly, sum);
    }

    #[test]
    fn external_product_trivial() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;

        for _ in 0..5 {
            let scalar = Poly::<TESTTYPE>::rand(1 << params.log_deg_ntru, rng);
            let x = NtruScalarCiphertext::trivial(scalar.clone(), params);

            let vector = Poly::<TESTTYPE>::rand(1 << params.log_deg_ntru, rng);
            let y = NtruVectorCiphertext::trivial(vector.clone(), params);

            let z = x.external_product(&y, params);
            let poly_prod = NtruScalarCiphertext::trivial(
                Poly(slow_negacyclic_mult(&scalar.0, &vector.0, TESTTYPE::ZERO)),
                params,
            );
            assert_eq!(z.ct, poly_prod.ct);
        }
    }

    #[test]
    fn external_product_noisy() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;
        let len = 1 << params.log_deg_ntru;
        let (key, _) = crate::key::NTRUKey::new(params, rng);

        for _ in 0..5 {
            for bit in 0..=1 {
                // accumulator is a ternary plaintext
                let mut scalar = Poly::<TESTTYPE>::new(len);
                scalar.iter_mut().for_each(|x| *x = sample_ternary(rng));
                let x = NtruScalarCiphertext::trivial(scalar.clone(), params);

                // vector ciphertexts are ternary monomials
                let mut vector1 = Poly::<TESTTYPE>::new(len);
                vector1.0[rng.gen_range(0..len)] = sample_ternary(rng);
                let y = NtruVectorCiphertext::trivial(vector1.clone(), params);

                let mut vector2 = Poly::<TESTTYPE>::new(len);
                vector2.0[rng.gen_range(0..len)] = sample_ternary(rng);
                let z = NtruVectorCiphertext::trivial(vector2.clone(), params);

                let mut bit_vec = vec![TESTTYPE::ZERO; len];
                bit_vec[0] = int_to_field(bit.into());
                let b = key.enc_bit_vec(bit, params, rng);

                let res = x
                    .external_product(&y, params)
                    .external_product(&z, params)
                    .external_product(&b, params);
                let dec = key.dec_scalar(res, params);
                let poly_prod = slow_negacyclic_mult(
                    &slow_negacyclic_mult(
                        &slow_negacyclic_mult(&scalar.0, &vector1.0, TESTTYPE::ZERO),
                        &vector2.0,
                        TESTTYPE::ZERO,
                    ),
                    &bit_vec,
                    TESTTYPE::ZERO,
                );
                assert_eq!(dec, Poly(poly_prod));
            }
        }
    }
}

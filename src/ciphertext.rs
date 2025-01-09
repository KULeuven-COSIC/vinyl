use swanky_field::PrimeFiniteField;

use crate::poly::{FFTPoly, Poly};

#[derive(Debug)]
pub struct LweCiphertext<M> {
    pub(crate) a: Vec<M>,
    pub(crate) b: M,
}

#[derive(Debug)]
pub struct NtruScalarCiphertext<T> {
    pub(crate) ct: Poly<T>,
}

// TODO: approximate decomposition
#[derive(Debug)]
pub struct NtruVectorCiphertext {
    pub(crate) ct: Vec<FFTPoly>,
}

impl<T: swanky_field::PrimeFiniteField> NtruScalarCiphertext<T> {
    // TODO: again assumes everything fits into a single limb
    /// Gadget decomposition into `n` parts. Returns an a vec such that
    /// `sum(result[i] * B^i) = self`
    fn gadget_decomp<BaI>(self, params: &crate::params::Params<T, BaI>) -> Vec<FFTPoly> {
        let mut res = Vec::with_capacity(params.dim_ngs);
        let mut as_ints: Vec<_> = self
            .ct
            .0
            .into_iter()
            .map(|x| x.into_int::<1>().as_words()[0])
            .collect();

        let base_int = params.gadget_base.into_int::<1>().as_words()[0];
        debug_assert_eq!(base_int.count_ones(), 1);
        debug_assert!(base_int.pow(params.dim_ngs as u32) > T::modulus_int::<1>().as_words()[0]);
        let base_mask = base_int - 1;
        let base_shift = base_int.ilog2();

        for _ in 0..params.dim_ngs {
            let mut part = Vec::<T>::with_capacity(as_ints.len());

            for coeff in as_ints.iter_mut() {
                part.push(
                    T::try_from_int::<1>((*coeff & base_mask).into())
                        .expect("Gadget decomp base guaranteed to fit in modulus"),
                );
                *coeff >>= base_shift;
            }

            res.push(params.fft.fwd(Poly(part)))
        }
        res
    }

    /// Compute the external product with a vector encryption
    /// Not implemented as a Mul trait for two reasons:
    /// - We need access to the params
    /// - This is costly, so explicit is better than implicit here
    pub(crate) fn external_product<BaI>(
        self,
        rhs: &NtruVectorCiphertext,
        params: &crate::params::Params<T, BaI>,
    ) -> Self {
        Self {
            ct: self
                .gadget_decomp(params)
                .into_iter()
                .zip(rhs.ct.iter())
                .map(|(x, y)| params.fft.inv(x * y))
                .sum(),
        }
    }
}

impl<T: swanky_field::PrimeFiniteField + for<'a> std::ops::Mul<&'a T, Output = T>>
    NtruScalarCiphertext<T>
{
    #[cfg(test)]
    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the "normal" scale $\Delta = Q/4$
    pub(crate) fn trivial<BaI>(pt: Poly<T>, params: &crate::params::Params<T, BaI>) -> Self {
        Self::trivial_scale(pt, &params.scale_ntru)
    }

    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the "halved" scale $\Delta = Q/8$ as used for the test vector
    pub(crate) fn trivial_half<BaI>(pt: Poly<T>, params: &crate::params::Params<T, BaI>) -> Self {
        Self::trivial_scale(pt, &params.half_scale_ntru)
    }

    /// Construct a trivial ciphertext embedding `pt`, with no noise;
    /// using the provided scale $\Delta$
    fn trivial_scale(pt: Poly<T>, scale: &T) -> Self {
        Self { ct: pt * scale }
    }
}

impl NtruVectorCiphertext {
    pub(crate) fn trivial<T: PrimeFiniteField, BaI>(
        m: Poly<T>,
        params: &crate::params::Params<T, BaI>,
    ) -> Self {
        let mut ct = Vec::with_capacity(params.dim_ngs);

        let mut base = T::ONE;
        for _ in 0..params.dim_ngs {
            ct.push(params.fft.fwd(m.clone() * base));
            base = base * params.gadget_base;
        }

        Self { ct }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::params::{TESTPARAMS, TESTTYPE};
    use crate::test_utils::*;

    use rand::Rng;
    use swanky_field::FiniteRing;

    #[test]
    fn gadget_decomp_consistency() {
        let rng = &mut rng(None);
        let poly = Poly::<TESTTYPE>::rand(1 << TESTPARAMS.log_deg_ntru, rng);
        let decomp = NtruScalarCiphertext { ct: poly.clone() }.gadget_decomp(&TESTPARAMS);
        let sum = decomp
            .into_iter()
            .enumerate()
            .map(|(i, x)| TESTPARAMS.fft.inv(x) * TESTPARAMS.gadget_base.pow_var_time(i as u128))
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
            let via_polys = NtruScalarCiphertext::trivial(
                Poly(slow_negacyclic_mult(&scalar.0, &vector.0, TESTTYPE::ZERO)),
                params,
            );
            assert_eq!(z.ct, via_polys.ct);
        }
    }

    #[test]
    fn external_product_noisy() {
        todo!()
    }
}

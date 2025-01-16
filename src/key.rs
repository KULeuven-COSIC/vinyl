//! Defining all sorts of cryptographic keys and the operations they can perform

use crate::ciphertext::{LweCiphertext, NtruScalarCiphertext, NtruVectorCiphertext};
use crate::modular::{
    bit_to_field, int_to_field, sample_discrete_gaussian, sample_ternary, ModSwitch, Modular,
};
use crate::params::{Params, Rng};
use crate::poly::{FFTPoly, Poly};

use crypto_bigint::subtle::ConditionallySelectable;
use swanky_field::{FiniteField, FiniteRing, PrimeFiniteField};

type LweKeyEl = u64;

/// An LWE secret key (binary)
#[derive(Debug)]
pub struct LWEKey {
    /// The actual key, stored in 64-bit chunks, from LSB to MSB
    key: Vec<LweKeyEl>,

    #[cfg(debug_assertions)]
    /// If we're running in debug mode, keep track of the LWE dimension for sanity checks
    dim: usize,
}

// TODO: We sometimes assume the modulus fits in a single limb
impl LWEKey {
    /// Sample a new random binary LWE key
    pub fn new(dim: usize, rng: &mut impl Rng) -> Self {
        // We read slightly more than we really have to
        // but we can ignore the remaining parts where needed
        let mut key = vec![0; dim.div_ceil(LweKeyEl::BITS as usize)];
        for x in key.iter_mut() {
            *x = rng.gen();
        }
        Self {
            key,
            #[cfg(debug_assertions)]
            dim,
        }
    }

    /// Get a bit, u8 because it's more of a number than a truth value
    fn get_bit(&self, i: usize) -> u8 {
        #[cfg(debug_assertions)]
        debug_assert!((0..self.dim).contains(&i));

        // Entirely little endian, to the bit level
        ((self.key[i / LweKeyEl::BITS as usize] >> (i % LweKeyEl::BITS as usize)) & 1) as u8
    }

    /// Obtain the `i`th component of the key
    /// as a bit (0 or 1 in the field).
    pub(crate) fn get<F: FiniteField>(&self, i: usize) -> F {
        bit_to_field(self.get_bit(i))
    }

    /// Iterate over the key bits
    /// Takes the LWE dimension as argument. This must match the params used to create the key!
    pub(crate) fn iter<F: FiniteField>(&self, dim: usize) -> impl Iterator<Item = F> + '_ {
        #[cfg(debug_assertions)]
        debug_assert_eq!(dim, self.dim);

        (0..dim).map(|i| self.get(i))
    }

    pub(crate) fn sample<F: PrimeFiniteField>(
        &self,
        dim: usize,
        stdev: f64,
        rng: &mut impl Rng,
    ) -> LweCiphertext<F> {
        let mut ct = LweCiphertext {
            a: vec![F::ZERO; dim],
            b: F::ZERO,
        };

        for (a, s) in ct.a.iter_mut().zip(self.iter(dim)) {
            *a = F::random(rng);
            ct.b -= *a * s;
        }
        ct.b += int_to_field(sample_discrete_gaussian(stdev, rng).into());

        ct
    }

    /// Encrypt a message
    pub fn encrypt<P: Params>(&self, message: u8, rng: &mut impl Rng) -> LweCiphertext<P::BaseInt>
    where
        P::BaseInt: PrimeFiniteField,
    {
        debug_assert!(message == 0 || message == 1);

        let mut ct = self.sample(P::DIM_LWE, P::ERR_STDEV_LWE, rng);
        let message =
            P::BaseInt::conditional_select(&P::BaseInt::ZERO, &P::BaseInt::ONE, message.into());
        ct.b += message * P::scale_lwe();

        ct
    }

    /// Decrypt a ciphertext with explicit parameters, rather than reading from a `Params`
    /// Useful e.g. when dealing with a different field
    fn decrypt_explicit<F: PrimeFiniteField>(
        &self,
        ct: LweCiphertext<F>,
        dim: usize,
        half_scale: F,
        scale: F,
    ) -> u8 {
        #[cfg(debug_assertions)]
        {
            debug_assert_eq!(ct.a.len(), dim);
            debug_assert_eq!(dim, self.dim);
            debug_assert!(self.key.len() * LweKeyEl::BITS as usize >= dim);
        }
        let mask: F = self.iter(dim).zip(ct.a).map(|(x, y): (F, F)| x * y).sum();

        // Add scale / 2 and floor div to round
        let out = (ct.b + mask + half_scale).into_int::<1>()
            / crypto_bigint::NonZero::new(scale.into_int::<1>()).unwrap();
        out.bit_vartime(0) as u8
    }

    /// Decrypt a ciphertext
    pub fn decrypt<P: Params>(&self, ct: LweCiphertext<P::BaseInt>) -> u8
    where
        P::BaseInt: PrimeFiniteField,
    {
        self.decrypt_explicit(ct, P::DIM_LWE, P::half_scale_lwe(), P::scale_lwe())
    }
}

/// An NTRU/NGS key, storing both the original polynomial 4f' + 1,
/// and its inverse mod X^N + 1
#[derive(Debug, Clone)]
pub struct NTRUKey {
    #[cfg(test)]
    /// `4*f' + 1`, the polynomial before inversion
    f: FFTPoly,
    /// `f^-1` the actual thing used for encryption in practice
    finv: FFTPoly,
}

impl NTRUKey {
    /// Generate a fresh NTRU/NGS key of given degree
    /// Also returns the coefficient vector that can be used to generate a NTRU->LWE KSK
    pub fn new<P: Params>(rng: &mut impl Rng) -> (Self, Poly<P::BootInt>)
    where
        P::BootInt: PrimeFiniteField,
    {
        let mut f = Poly::<P::BootInt>::new(1 << P::LOG_DEG_NTRU);

        // Until we find an invertible vector
        loop {
            f.iter_mut()
                .for_each(|x| *x = P::scale_ntru_key() * sample_ternary(rng));
            f = f + P::BootInt::ONE;
            if let Some(finv) = f.invert() {
                return (
                    Self {
                        #[cfg(test)]
                        f: P::fft().fwd(f.clone()),
                        finv: P::fft().fwd(finv),
                    },
                    f,
                );
            }
        }
    }

    pub fn enc_bit_vec<P: Params>(&self, bit: u8, rng: &mut impl Rng) -> NtruVectorCiphertext
    where
        P::BootInt: PrimeFiniteField,
    {
        let fft = P::fft();
        let mut res = Vec::with_capacity(P::DIM_NGS);
        let mut gadget_base_pow =
            P::BootInt::try_from_int::<1>(bit.into()).expect("A bit should fit any field");
        for _ in 0..P::DIM_NGS {
            let mut g = Poly::<P::BootInt>::new(1 << P::LOG_DEG_NTRU);
            g.iter_mut().for_each(|x| {
                // TODO: is this discrete gaussian or ternary?
                *x = int_to_field(sample_discrete_gaussian(P::ERR_STDEV_NTRU, rng).into())
            });

            // g / f
            let mut component: Poly<P::BootInt> = fft.inv::<P::BootInt>(fft.fwd(g) * &self.finv);
            // g / f + m * B^i
            component.0[0] += gadget_base_pow;

            gadget_base_pow *= P::gadget_base();
            res.push(fft.fwd(component))
        }
        NtruVectorCiphertext { ct: res }
    }

    #[cfg(test)]
    pub(crate) fn dec_scalar<P: Params>(
        &self,
        ct: NtruScalarCiphertext<P::BootInt>,
    ) -> Poly<P::BootInt>
    where
        P::BootInt: PrimeFiniteField,
    {
        // f * (g/f + Δ m) = g + f Δ m = g' + Q/Δ f' Δ  m + Δ m = g' + Δ m

        use crate::modular::lift_centered;
        let delta_m_plus_noise = P::fft().inv::<P::BootInt>(P::fft().fwd(ct.ct) * &self.f);
        let modulus = P::BootInt::modulus_int::<1>().as_words()[0] as i64;

        Poly(
            delta_m_plus_noise
                .0
                .into_iter()
                .map(|coef| {
                    let mut ccoef = lift_centered(coef);
                    if ccoef < 0 {
                        ccoef -= lift_centered(P::half_scale_ntru());
                    } else {
                        ccoef += lift_centered(P::half_scale_ntru());
                    }
                    let intval = ccoef / lift_centered(P::scale_ntru());
                    P::BootInt::try_from_int::<1>((((intval + modulus) % modulus) as u64).into())
                        .unwrap()
                })
                .collect(),
        )
    }
}

// TODO: Generalize to mkLWE
/// A key-switching key from NTRU to LWE
// Outer goes over the decomposition, inner over coefficients of f
#[derive(Clone, Debug)]
struct KskNtruLwe<F>(Vec<Vec<LweCiphertext<F>>>);

impl<F: PrimeFiniteField> KskNtruLwe<F> {
    /// Build a KSK that transfers from ciphertext under the ntru key
    /// to ciphertexts under the lwe key
    fn new<P: Params<BootInt = F>>(ntru: &Poly<F>, lwe: &LWEKey, rng: &mut impl Rng) -> Self {
        let mut sample = |delta: F| -> LweCiphertext<F> {
            let mut ct = lwe.sample(P::DIM_LWE, P::ERR_STDEV_LWE, rng);
            debug_assert_ne!(ct.a.len(), 0);
            ct.b += delta;
            ct
        };

        let mut res = Vec::with_capacity(P::KSK_NTRU_LWE_DIM);
        let mut base = F::ONE;
        for _ in 0..P::KSK_NTRU_LWE_DIM {
            let mut row = Vec::with_capacity(ntru.degree() + 1);
            row.push(sample(base * ntru.0[0]));
            for coef in ntru.0.iter().skip(1).rev() {
                row.push(sample(-base * *coef));
            }
            res.push(row);
            base *= P::ksk_ntru_lwe_base();
        }
        Self(res)
    }

    pub fn key_switch<P: Params<BootInt = F>>(
        &self,
        ct: NtruScalarCiphertext<F>,
    ) -> LweCiphertext<F> {
        let decomp = ct.gadget_decomp(P::KSK_NTRU_LWE_DIM, P::ksk_ntru_lwe_base());
        // Ugly roundtrip through Poly so that we get vectorial addition for free
        let mut a = Poly::new(P::DIM_LWE);
        let mut b = F::ZERO;

        // First on the power of the gadget decomp
        #[cfg(debug_assertions)]
        debug_assert_eq!(decomp.len(), self.0.len());
        for (y, a_decomp) in decomp.into_iter().zip(self.0.iter()) {
            #[cfg(debug_assertions)]
            debug_assert_eq!(y.0.len(), a_decomp.len());
            for (coef, sample) in y.0.into_iter().zip(a_decomp) {
                a = a + Poly(sample.a.clone()) * coef;
                b += sample.b * coef;
            }
        }
        LweCiphertext { a: a.0, b }
    }
}

/// A FINAL bootstrapping key: encryptions of an LWE secret key under the NTRU boot key
/// ready for use in CMUX gates.
#[derive(Clone, Debug)]
struct Bsk(Vec<NtruVectorCiphertext>);

impl Bsk {
    fn new<P: Params>(boot_key: &NTRUKey, base_key: &LWEKey, rng: &mut impl Rng) -> Self
    where
        P::BootInt: PrimeFiniteField,
    {
        Self(Vec::from_iter((0..P::DIM_LWE).map(|idx| {
            boot_key.enc_bit_vec::<P>(base_key.get_bit(idx), rng)
        })))
    }

    /// Compute `acc` * CMUX_`i`(`a`) = `acc` * VEnc(1 + (X^(`a`) - 1) * bsk_`i`)
    fn ext_prod_cmux<P: Params>(
        &self,
        acc: NtruScalarCiphertext<P::BootInt>,
        i: usize,
        a: P::ExpRing,
        fft: &crate::fft::FFTPlan,
    ) -> NtruScalarCiphertext<P::BootInt>
    where
        P::BootInt: PrimeFiniteField,
        <P::ExpRing as Modular>::Single: TryInto<usize>,
        <<P::ExpRing as Modular>::Single as TryInto<usize>>::Error: std::fmt::Debug,
    {
        // acc * (1 + (X^a - 1) * bsk) = acc + (acc * X^a - acc) * bsk
        // = acc + (rot(acc, a) - acc) * bsk
        let orig = acc.clone();
        orig.clone()
            + (acc.rot(
                a.extract()
                    .try_into()
                    .expect("coefficient doesn't fit usize"),
            ) - orig)
                .external_product::<P>(&self.0[i], fft)
    }
}

/// A FINAL public key, contains keyswitching and bootstrapping keys
#[derive(Clone, Debug)]
pub struct PublicKey<P: Params> {
    ksk: KskNtruLwe<P::BootInt>,
    bsk: Bsk,
}

impl<P: Params> PublicKey<P>
where
    P::BaseInt: ModSwitch<P::ExpRing>,
    P::BootInt: ModSwitch<P::BaseInt>
        + PrimeFiniteField
        + for<'a> std::ops::Mul<&'a P::BootInt, Output = P::BootInt>,
    <P::ExpRing as Modular>::Single: TryInto<usize>,
    <<P::ExpRing as Modular>::Single as TryInto<usize>>::Error: std::fmt::Debug,
{
    pub fn bootstrap(&self, ct: LweCiphertext<P::BaseInt>) -> LweCiphertext<P::BaseInt> {
        let fft = P::fft();
        let ct_er = ct.modswitch();
        let mut test_vector = Poly::new(1 << P::LOG_DEG_NTRU);
        test_vector.0.iter_mut().for_each(|x| *x = P::BootInt::ONE);
        let mut acc = NtruScalarCiphertext::trivial_half::<P>(test_vector).rot(
            ct_er
                .b
                .extract()
                .try_into()
                .expect("Weird platform bit sizes??")
                + (1 << (P::LOG_DEG_NTRU - 1)),
        );

        for (i, a) in ct_er.a.into_iter().enumerate() {
            acc = self.bsk.ext_prod_cmux::<P>(acc, i, a, fft);
        }

        // add Q/8 * sum_i X^i
        acc.ct.iter_mut().for_each(|x| *x += P::half_scale_ntru());

        let acc_lwe = self.ksk.key_switch::<P>(acc);
        acc_lwe.modswitch()
    }
}

/// A full FINAL key, including all private key material
#[derive(Debug)]
pub struct Key<P: Params> {
    base: LWEKey,
    #[cfg(test)]
    /// The key used for the NGS things, decryption isn't needed during regular execution
    /// but to debug could be useful
    boot: NTRUKey,
    public: PublicKey<P>,
}

impl<P: Params> Key<P>
where
    P::BootInt: PrimeFiniteField,
    P::BaseInt: PrimeFiniteField,
{
    pub fn new(rng: &mut impl Rng) -> Self {
        let base = LWEKey::new(P::DIM_LWE, rng);
        let (boot, boot_coefs) = NTRUKey::new::<P>(rng);
        let ksk = KskNtruLwe::new::<P>(&boot_coefs, &base, rng);
        let bsk = Bsk::new::<P>(&boot, &base, rng);

        Key {
            base,
            #[cfg(test)]
            boot,
            public: PublicKey { ksk, bsk },
        }
    }

    pub fn export(&self) -> &PublicKey<P> {
        &self.public
    }

    pub fn encrypt(&self, val: u8, rng: &mut impl Rng) -> LweCiphertext<P::BaseInt> {
        self.base.encrypt::<P>(val, rng)
    }

    pub fn decrypt(&self, ct: LweCiphertext<P::BaseInt>) -> u8 {
        self.base.decrypt::<P>(ct)
    }
}

#[cfg(test)]
mod test {
    use rand::Rng;

    use super::*;
    use crate::{params::TestParams, test_utils::*};
    type BootInt = <TestParams as Params>::BootInt;

    #[test]
    fn lwe_roundtrip() {
        let rng = &mut rng(None);
        let key = LWEKey::new(TestParams::DIM_LWE, rng);
        for _ in 0..100 {
            let ct0 = key.encrypt::<TestParams>(0, rng);
            assert_eq!(key.decrypt::<TestParams>(ct0), 0);
            let ct1 = key.encrypt::<TestParams>(1, rng);
            assert_eq!(key.decrypt::<TestParams>(ct1), 1);
        }
    }

    #[test]
    fn ntru_decrypt_trivial() {
        let rng = &mut rng(None);
        let len = 1 << TestParams::LOG_DEG_NTRU;
        let (key, _) = crate::key::NTRUKey::new::<TestParams>(rng);

        let mut scalar = Poly::<BootInt>::new(len);
        scalar.iter_mut().for_each(|x| *x = sample_ternary(rng));
        let x = NtruScalarCiphertext::trivial::<TestParams>(scalar.clone());

        assert_eq!(scalar, key.dec_scalar::<TestParams>(x));
    }

    #[test]
    fn ntru_decrypt_noisy() {
        let rng = &mut rng(None);
        let len = 1 << TestParams::LOG_DEG_NTRU;

        // We'll try a few random bit encryptions, but those are vector ciphertexts
        // so we need this scalar to bring it back into scalar with an external product
        let zero = Poly::<BootInt>::new(len);
        let one = Poly::new(len) + BootInt::ONE;
        let acc = NtruScalarCiphertext::trivial::<TestParams>(one.clone());
        let fft = TestParams::fft();

        for _ in 0..5 {
            let (key, _) = NTRUKey::new::<TestParams>(rng);
            let ct = key.enc_bit_vec::<TestParams>(0, rng);
            assert_eq!(
                key.dec_scalar::<TestParams>(acc.clone().external_product::<TestParams>(&ct, fft)),
                zero
            );
            let ct = key.enc_bit_vec::<TestParams>(1, rng);
            assert_eq!(
                key.dec_scalar::<TestParams>(acc.clone().external_product::<TestParams>(&ct, fft)),
                one
            );
        }
    }

    #[test]
    fn trivial_ksk_dec() {
        let rng = &mut rng(None);
        let lwe_key = LWEKey::new(TestParams::DIM_LWE, rng);
        let (_, ntru_coefs) = NTRUKey::new::<TestParams>(rng);
        let ksk = KskNtruLwe::new::<TestParams>(&ntru_coefs, &lwe_key, rng);

        for _ in 0..100 {
            let mut plain = Poly::new(1 << TestParams::LOG_DEG_NTRU);
            plain
                .iter_mut()
                .for_each(|x| *x = int_to_field(rng.gen_range(0..=1u8).into()));

            let ntru_ct = NtruScalarCiphertext::trivial::<TestParams>(plain.clone());
            let lwe_ct = ksk.key_switch::<TestParams>(ntru_ct);

            let dec = lwe_key.decrypt_explicit(
                lwe_ct,
                TestParams::DIM_LWE,
                TestParams::half_scale_ntru(),
                TestParams::scale_ntru(),
            );
            assert_eq!(plain.0[0], int_to_field(dec.into()));
        }
    }

    #[test]
    fn enc_ksk_dec() {
        let rng = &mut rng(None);
        let lwe_key = LWEKey::new(TestParams::DIM_LWE, rng);
        let (ntru_key, ntru_coefs) = NTRUKey::new::<TestParams>(rng);
        let scalar_one = NtruScalarCiphertext::trivial::<TestParams>(
            Poly::new(1 << TestParams::LOG_DEG_NTRU) + BootInt::ONE,
        );
        let ksk = KskNtruLwe::new::<TestParams>(&ntru_coefs, &lwe_key, rng);
        let fft = TestParams::fft();

        for _ in 0..100 {
            for b in 0..=1 {
                let ntru_ct = scalar_one.clone().external_product::<TestParams>(
                    &ntru_key.enc_bit_vec::<TestParams>(b, rng),
                    fft,
                );

                assert_eq!(
                    ntru_key.dec_scalar::<TestParams>(ntru_ct.clone()).0[0],
                    int_to_field(b.into())
                );

                let lwe_ct = ksk.key_switch::<TestParams>(ntru_ct);
                let dec = lwe_key.decrypt_explicit(
                    lwe_ct,
                    TestParams::DIM_LWE,
                    TestParams::half_scale_ntru(),
                    TestParams::scale_ntru(),
                );
                assert_eq!(b, dec);
            }
        }
    }

    #[test]
    fn bootstrap() {
        let rng = &mut rng(None);
        let key = Key::<TestParams>::new(rng);
        let evk = key.export();

        for _ in 0..100 {
            for b in 0..=1 {
                let ct = key.base.encrypt::<TestParams>(b, rng);
                assert_eq!(key.base.decrypt::<TestParams>(ct.clone()), b);
                let ct1 = evk.bootstrap(ct);
                assert_eq!(key.base.decrypt::<TestParams>(ct1.clone()), b);
                let ct2 = evk.bootstrap(ct1);
                assert_eq!(key.base.decrypt::<TestParams>(ct2), b);
            }
        }
    }
}

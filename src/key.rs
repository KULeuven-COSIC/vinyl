//! Defining all sorts of cryptographic keys and the operations they can perform

use crate::ciphertext::{LweCiphertext, NtruScalarCiphertext, NtruVectorCiphertext};
use crate::modular::{
    bit_to_field, int_to_field, sample_discrete_gaussian, sample_ternary, ModSwitch, Modular,
};
use crate::params::{Params, Rng};
use crate::poly::{FFTPoly, Poly};

use swanky_field::{FiniteField, PrimeFiniteField};

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
        let mut key = vec![0; (dim + LweKeyEl::BITS as usize - 1) / LweKeyEl::BITS as usize];
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
            ct.b += *a * s;
        }
        ct.b += int_to_field(sample_discrete_gaussian(stdev, rng).into());

        ct
    }

    /// Encrypt a message
    pub fn encrypt<BoI, BaI, ER>(
        &self,
        message: u8,
        params: &Params<BoI, BaI, ER>,
        rng: &mut impl Rng,
    ) -> LweCiphertext<BaI>
    where
        BaI: PrimeFiniteField,
    {
        debug_assert!(message == 0 || message == 1);

        let mut ct = self.sample(params.dim_lwe, params.err_stdev_lwe, rng);
        let message = BaI::conditional_select(&BaI::ZERO, &BaI::ONE, message.into());
        ct.b += message * params.scale_lwe;

        ct
    }

    /// Decrypt a ciphertext with explicit parameters, rather than reading from a `Params`
    /// Useful e.g. when dealing with a different field
    fn decrypt_explicit<F: PrimeFiniteField + for<'a> std::ops::Add<&'a F, Output = F>>(
        &self,
        ct: LweCiphertext<F>,
        dim: usize,
        half_scale: &F,
        scale: &F,
    ) -> u8 {
        #[cfg(debug_assertions)]
        {
            debug_assert_eq!(ct.a.len(), dim);
            debug_assert_eq!(dim, self.dim);
            debug_assert!(self.key.len() * LweKeyEl::BITS as usize >= dim);
        }
        let mask = self.iter(dim).zip(ct.a).map(|(x, y): (F, F)| x * y).sum();

        // Add scale / 2 and floor div to round
        let out = (ct.b - mask + half_scale).into_int::<1>()
            / crypto_bigint::NonZero::new(scale.into_int::<1>()).unwrap();
        out.bit_vartime(0) as u8
    }

    /// Decrypt a ciphertext
    pub fn decrypt<BoI, BaI, ER>(&self, ct: LweCiphertext<BaI>, params: &Params<BoI, BaI, ER>) -> u8
    where
        BaI: PrimeFiniteField + for<'a> std::ops::Add<&'a BaI, Output = BaI>,
    {
        self.decrypt_explicit(
            ct,
            params.dim_lwe,
            &params.half_scale_lwe,
            &params.scale_lwe,
        )
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
    pub fn new<BoI: PrimeFiniteField, BaI, ER>(
        params: &Params<BoI, BaI, ER>,
        rng: &mut impl Rng,
    ) -> (Self, Poly<BoI>) {
        let mut f = Poly::<BoI>::new(1 << params.log_deg_ntru);

        // Until we find an invertible vector
        loop {
            f.iter_mut()
                .for_each(|x| *x = params.scale_ntru_key * sample_ternary(rng));
            f = f + BoI::ONE;
            if let Some(finv) = f.invert() {
                return (
                    Self {
                        #[cfg(test)]
                        f: params.fft.fwd(f.clone()),
                        finv: params.fft.fwd(finv),
                    },
                    f,
                );
            }
        }
    }

    pub fn enc_bit_vec<BoI: PrimeFiniteField, BaI, ER>(
        &self,
        bit: u8,
        params: &Params<BoI, BaI, ER>,
        rng: &mut impl Rng,
    ) -> NtruVectorCiphertext {
        let mut res = Vec::with_capacity(params.dim_ngs);
        let mut gadget_base_pow =
            BoI::try_from_int::<1>(bit.into()).expect("A bit should fit any field");
        for _ in 0..params.dim_ngs {
            let mut g = Poly::<BoI>::new(1 << params.log_deg_ntru);
            g.iter_mut().for_each(|x| {
                // TODO: is this discrete gaussian or ternary?
                *x = int_to_field(sample_discrete_gaussian(params.err_stdev_ntru, rng).into())
            });

            // g / f
            let mut component: Poly<BoI> = params.fft.inv::<BoI>(params.fft.fwd(g) * &self.finv);
            // g / f + m * B^i
            component.0[0] += gadget_base_pow;

            gadget_base_pow *= params.gadget_base;
            res.push(params.fft.fwd(component))
        }
        NtruVectorCiphertext { ct: res }
    }

    #[cfg(test)]
    pub(crate) fn dec_scalar<BoI: PrimeFiniteField, BaI, ER>(
        &self,
        ct: NtruScalarCiphertext<BoI>,
        params: &Params<BoI, BaI, ER>,
    ) -> Poly<BoI> {
        // f * (g/f + Δ m) = g + f Δ m = g' + Q/Δ f' Δ  m + Δ m = g' + Δ m

        use crate::modular::lift_centered;
        let delta_m_plus_noise = params.fft.inv::<BoI>(params.fft.fwd(ct.ct) * &self.f);
        let modulus = BoI::modulus_int::<1>().as_words()[0] as i64;

        Poly(
            delta_m_plus_noise
                .0
                .into_iter()
                .map(|coef| {
                    let mut ccoef = lift_centered(coef);
                    if ccoef < 0 {
                        ccoef -= lift_centered(params.half_scale_ntru);
                    } else {
                        ccoef += lift_centered(params.half_scale_ntru);
                    }
                    let intval = ccoef / lift_centered(params.scale_ntru);
                    BoI::try_from_int::<1>((((intval + modulus) % modulus) as u64).into()).unwrap()
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
    fn new<BaI, ER>(
        ntru: &Poly<F>,
        lwe: &LWEKey,
        params: &Params<F, BaI, ER>,
        rng: &mut impl Rng,
    ) -> Self {
        let mut res = Vec::with_capacity(params.ksk_ntru_lwe_dim);
        let mut base = F::ONE;
        for _ in 0..params.ksk_ntru_lwe_dim {
            let mut row = Vec::with_capacity(ntru.degree() + 1);
            for coef in ntru.0.iter() {
                let mut ct = lwe.sample(params.dim_lwe, params.err_stdev_lwe, rng);
                debug_assert_ne!(ct.a.len(), 0);
                ct.b += base.clone() * *coef;
                row.push(ct)
            }
            res.push(row);
            base = base * params.ksk_ntru_lwe_base;
        }
        Self(res)
    }

    pub fn key_switch<BaI, ER>(
        &self,
        ct: NtruScalarCiphertext<F>,
        params: &Params<F, BaI, ER>,
    ) -> LweCiphertext<F> {
        let decomp = ct.gadget_decomp(params.ksk_ntru_lwe_dim, params.ksk_ntru_lwe_base);
        // Ugly roundtrip through Poly so that we get vectorial addition for free
        let mut a = Poly::new(params.dim_lwe);
        let mut b = F::ZERO;

        // First on the power of the gadget decomp
        #[cfg(debug_assertions)]
        debug_assert_eq!(decomp.len(), self.0.len());
        for (y, a_decomp) in decomp.into_iter().zip(self.0.iter()) {
            #[cfg(debug_assertions)]
            debug_assert_eq!(y.0.len(), a_decomp.len());
            for (coef, sample) in y.0.into_iter().zip(a_decomp) {
                a = a + Poly(sample.a.clone()) * coef;
                b = b + sample.b * coef;
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
    fn new<BoI: PrimeFiniteField, BaI, ER>(
        boot_key: &NTRUKey,
        base_key: &LWEKey,
        params: &Params<BoI, BaI, ER>,
        rng: &mut impl Rng,
    ) -> Self {
        Self(Vec::from_iter((0..params.dim_lwe).map(|idx| {
            boot_key.enc_bit_vec(base_key.get_bit(idx), params, rng)
        })))
    }

    /// Compute `acc` * CMUX_`i`(`a`) = `acc` * VEnc(1 + (X^(`a`) - 1) * bsk_`i`)
    fn ext_prod_cmux<BoI, BaI, ER>(
        &self,
        acc: NtruScalarCiphertext<BoI>,
        i: usize,
        a: ER,
        params: &crate::params::Params<BoI, BaI, ER>,
    ) -> NtruScalarCiphertext<BoI>
    where
        BoI: PrimeFiniteField,
        ER: Modular,
        ER::Single: TryInto<usize>,
        <ER::Single as TryInto<usize>>::Error: std::fmt::Debug,
    {
        // acc * (1 + (X^a - 1) * bsk) = acc + (acc * X^a - acc) * bsk
        // = acc + (rot(acc, a) - acc) * bsk
        let orig = acc.clone();
        (acc.rot(
            a.extract()
                .try_into()
                .expect("coefficient doesn't fit usize"),
        ) - orig)
            .external_product(&self.0[i], params)
    }
}

/// A FINAL public key, contains keyswitching and bootstrapping keys
#[derive(Clone, Debug)]
pub struct PublicKey<'a, BootInt, BaseInt, ER> {
    params: &'a Params<BootInt, BaseInt, ER>,
    ksk: KskNtruLwe<BootInt>,
    bsk: Bsk,
}

impl<'a, BoI, BaI, ER> PublicKey<'a, BoI, BaI, ER>
where
    BaI: ModSwitch<ER>,
    BoI: ModSwitch<BaI> + PrimeFiniteField + for<'b> std::ops::Mul<&'b BoI, Output = BoI>,
    ER: Modular,
    ER::Single: TryInto<usize>,
    <ER::Single as TryInto<usize>>::Error: std::fmt::Debug,
{
    pub fn bootstrap(
        &self,
        ct: LweCiphertext<BaI>,
        params: &Params<BoI, BaI, ER>,
    ) -> LweCiphertext<BaI> {
        let ct_er = ct.modswitch();
        let mut test_vector = Poly::new(1 << params.log_deg_ntru);
        test_vector.0.iter_mut().for_each(|x| *x = BoI::ONE);
        let mut acc = NtruScalarCiphertext::trivial_half(test_vector, params).external_product(
            &NtruVectorCiphertext::monomial(
                ct_er
                    .b
                    .extract()
                    .try_into()
                    .expect("Weird platform bit sizes??")
                    + (1 << (params.log_deg_ntru - 1)),
                params,
            ),
            params,
        );

        for (i, a) in ct_er.a.into_iter().enumerate() {
            acc = self.bsk.ext_prod_cmux(acc, i, a, params);
        }

        // add Q/8 * sum_i X^i
        acc.ct
            .iter_mut()
            .for_each(|x| *x = *x + params.half_scale_ntru);

        let acc_lwe = self.ksk.key_switch(acc, params);
        acc_lwe.modswitch()
    }
}

/// A full FINAL key, including all private key material
#[derive(Debug)]
pub struct Key<'a, BootInt, BaseInt, ER> {
    params: &'a Params<BootInt, BaseInt, ER>,
    base: LWEKey,
    #[cfg(debug_assertions)]
    /// The key used for the NGS things, decryption isn't needed during regular execution
    /// but to debug could be useful
    boot: NTRUKey,
    public: PublicKey<'a, BootInt, BaseInt, ER>,
}

impl<'a, BoI, BaI, ER> Key<'a, BoI, BaI, ER>
where
    BoI: PrimeFiniteField,
    BaI: PrimeFiniteField,
{
    pub fn new(params: &'a Params<BoI, BaI, ER>, rng: &mut impl Rng) -> Self {
        let base = LWEKey::new(params.dim_lwe, rng);
        let (boot, boot_coefs) = NTRUKey::new(params, rng);
        let ksk = KskNtruLwe::new(&boot_coefs, &base, params, rng);
        let bsk = Bsk::new(&boot, &base, params, rng);

        Key {
            params,
            base,
            #[cfg(debug_assertions)]
            boot,
            public: PublicKey { params, ksk, bsk },
        }
    }

    pub fn export(&self) -> &PublicKey<BoI, BaI, ER> {
        &self.public
    }
}

#[cfg(test)]
mod test {
    use rand::Rng;
    use swanky_field::FiniteRing;

    use super::*;
    use crate::{
        params::{TESTPARAMS, TESTTYPE},
        test_utils::*,
    };

    #[test]
    fn lwe_roundtrip() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;
        let key = LWEKey::new(params.dim_lwe, rng);
        for _ in 0..100 {
            let ct0 = key.encrypt(0, params, rng);
            assert_eq!(key.decrypt(ct0, params), 0);
            let ct1 = key.encrypt(1, params, rng);
            assert_eq!(key.decrypt(ct1, params), 1);
        }
    }

    #[test]
    fn ntru_decrypt_trivial() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;
        let len = 1 << params.log_deg_ntru;
        let (key, _) = crate::key::NTRUKey::new(params, rng);

        let mut scalar = Poly::<TESTTYPE>::new(len);
        scalar.iter_mut().for_each(|x| *x = sample_ternary(rng));
        let x = NtruScalarCiphertext::trivial(scalar.clone(), params);

        assert_eq!(scalar, key.dec_scalar(x, params));
    }

    #[test]
    fn ntru_decrypt_noisy() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;
        let len = 1 << params.log_deg_ntru;

        // We'll try a few random bit encryptions, but those are vector ciphertexts
        // so we need this scalar to bring it back into scalar with an external product
        let zero = Poly::new(len);
        let one = Poly::new(len) + TESTTYPE::ONE;
        let acc = NtruScalarCiphertext::trivial(one.clone(), params);

        for _ in 0..5 {
            let (key, _) = NTRUKey::new(params, rng);
            let ct = key.enc_bit_vec(0, params, rng);
            assert_eq!(
                key.dec_scalar(acc.clone().external_product(&ct, params), params),
                zero
            );
            let ct = key.enc_bit_vec(1, params, rng);
            assert_eq!(
                key.dec_scalar(acc.clone().external_product(&ct, params), params),
                one
            );
        }
    }

    #[test]
    fn enc_ksk_dec() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;
        let lwe_key = LWEKey::new(params.dim_lwe, rng);
        let (_, ntru_coefs) = NTRUKey::new(params, rng);
        let ksk = KskNtruLwe::new(&ntru_coefs, &lwe_key, params, rng);

        for _ in 0..100 {
            let mut plain = Poly::new(1 << params.log_deg_ntru);
            plain
                .iter_mut()
                .for_each(|x| *x = int_to_field(rng.gen_range(0..=1u8).into()));

            let ntru_ct = NtruScalarCiphertext::trivial(plain.clone(), params);
            let lwe_ct = ksk.key_switch(ntru_ct, params);

            let dec = lwe_key.decrypt_explicit(
                lwe_ct,
                params.dim_lwe,
                &params.half_scale_ntru,
                &params.scale_ntru,
            );
            assert_eq!(plain.0[0], int_to_field(dec.into()));
        }
    }

    #[test]
    fn bootstrap() {
        let rng = &mut rng(None);
        let params = &TESTPARAMS;
        let key = Key::new(&params, rng);
        let evk = key.export();

        for _ in 0..10 {
            for b in 0..=1 {
                let ct = key.base.encrypt(b, params, rng);
                let ct1 = evk.bootstrap(ct.clone(), params);
                let ct2 = evk.bootstrap(ct1.clone(), params);
                assert_eq!(key.base.decrypt(ct, params), b);
                assert_eq!(key.base.decrypt(ct1, params), b);
                assert_eq!(key.base.decrypt(ct2, params), b);
            }
        }
    }
}

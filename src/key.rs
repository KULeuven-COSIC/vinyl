//! Defining all sorts of cryptographic keys and the operations they can perform

use crate::ciphertext::{LweCiphertext, NtruScalarCiphertext, NtruVectorCiphertext};
use crate::modular::{bit_to_field, int_to_field, sample_discrete_gaussian, sample_ternary};
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
    pub fn new<BoI, BaI>(params: &Params<BoI, BaI>, rng: &mut impl Rng) -> Self {
        // We read slightly more than we really have to
        // but we can ignore the remaining parts where needed
        let mut key =
            vec![0; (params.dim_lwe + LweKeyEl::BITS as usize - 1) / LweKeyEl::BITS as usize];
        for x in key.iter_mut() {
            *x = rng.gen();
        }
        Self {
            key,
            #[cfg(debug_assertions)]
            dim: params.dim_lwe,
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

    /// Encrypt a message
    pub fn encrypt<BaI, BoI>(
        &self,
        message: u8,
        params: &Params<BoI, BaI>,
        rng: &mut impl Rng,
    ) -> LweCiphertext<BaI>
    where
        BaI: PrimeFiniteField,
    {
        debug_assert!(message == 0 || message == 1);

        let mut ct = LweCiphertext {
            a: Vec::with_capacity(params.dim_lwe),
            b: BaI::ZERO,
        };

        for (a, s) in ct.a.iter_mut().zip(self.iter(params.dim_lwe)) {
            *a = BaI::random(rng);
            ct.b += *a * s;
        }

        let message = BaI::conditional_select(&BaI::ZERO, &BaI::ONE, message.into());
        ct.b += message * params.scale_lwe
            + int_to_field(sample_discrete_gaussian(params.err_stdev_lwe, rng).into());

        ct
    }

    /// Decrypt a ciphertext
    pub fn decrypt<BaI, BoI>(&self, ct: LweCiphertext<BaI>, params: &Params<BoI, BaI>) -> u8
    where
        BaI: PrimeFiniteField,
    {
        let mask = self
            .iter(params.dim_lwe)
            .zip(ct.a)
            .map(|(x, y): (BaI, BaI)| x * y)
            .sum();

        // Add scale / 2 and floor div to round
        let out = (ct.b - mask + params.half_scale_lwe).into_int::<1>()
            / crypto_bigint::NonZero::new(params.scale_lwe.into_int::<1>()).unwrap();
        out.bit_vartime(0) as u8
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
    pub fn new<BoI: PrimeFiniteField, BaI>(params: &Params<BoI, BaI>, rng: &mut impl Rng) -> Self {
        let mut f = Poly::<BoI>::new(1 << params.log_deg_ntru);
        let scale: BoI = int_to_field(4u8.into());

        // Until we find an invertible vector
        loop {
            f.iter_mut().for_each(|x| *x = scale * sample_ternary(rng));
            f = f + BoI::ONE;
            if let Some(finv) = f.invert() {
                return Self {
                    #[cfg(test)]
                    f: params.fft.fwd(f),
                    finv: params.fft.fwd(finv),
                };
            }
        }
    }

    pub fn enc_bit_vec<BoI: PrimeFiniteField, BaI>(
        &self,
        bit: u8,
        params: &Params<BoI, BaI>,
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
    pub(crate) fn dec_scalar<BoI: PrimeFiniteField, BaI>(
        &self,
        ct: NtruScalarCiphertext<BoI>,
        params: &Params<BoI, BaI>,
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

// TODO
// also TODO: Distinguish from LWE->LWE KSK
/// A key-switching key from NTRU to LWE
#[derive(Clone, Debug)]
struct KskNtruLwe<F>(std::marker::PhantomData<F>);

impl<F> KskNtruLwe<F> {
    /// Build a KSK that transfers from ciphertext under the ntru key
    /// to ciphertexts under the lwe key
    fn new(ntru: &NTRUKey, lwe: &LWEKey) -> Self {
        todo!()
    }
}

/// A FINAL bootstrapping key: encryptions of an LWE secret key under the NTRU boot key
/// ready for use in CMUX gates.
#[derive(Clone, Debug)]
struct Bsk(Vec<NtruVectorCiphertext>);

impl Bsk {
    fn new<BoI: PrimeFiniteField, BaI>(
        boot_key: &NTRUKey,
        base_key: &LWEKey,
        params: &Params<BoI, BaI>,
        rng: &mut impl Rng,
    ) -> Self {
        Self(Vec::from_iter((0..params.dim_lwe).map(|idx| {
            boot_key.enc_bit_vec(base_key.get_bit(idx), params, rng)
        })))
    }
}

/// A FINAL public key, contains keyswitching and bootstrapping keys
#[derive(Clone, Debug)]
pub struct PublicKey<'a, BootInt, BaseInt> {
    params: &'a Params<BootInt, BaseInt>,
    ksk: KskNtruLwe<BootInt>,
    bsk: Bsk,
}

/// A full FINAL key, including all private key material
#[derive(Debug)]
pub struct Key<'a, BootInt, BaseInt> {
    params: &'a Params<BootInt, BaseInt>,
    base: LWEKey,
    #[cfg(debug_assertions)]
    /// The key used for the NGS things, decryption isn't needed during regular execution
    /// but to debug could be useful
    boot: NTRUKey,
    public: PublicKey<'a, BootInt, BaseInt>,
}

impl<'a, BoI, BaI> Key<'a, BoI, BaI>
where
    BoI: PrimeFiniteField,
    BaI: PrimeFiniteField,
{
    pub fn new(params: &'a Params<BoI, BaI>, rng: &mut impl Rng) -> Self {
        let base = LWEKey::new(params, rng);
        let boot = NTRUKey::new(params, rng);
        let ksk = KskNtruLwe::new(&boot.clone(), &base); // Check what kind of modswitching we may need
        let bsk = Bsk::new(&boot, &base, params, rng);

        Key {
            params,
            base,
            #[cfg(debug_assertions)]
            boot,
            public: PublicKey { params, ksk, bsk },
        }
    }

    pub fn export(&self) -> &PublicKey<BoI, BaI> {
        &self.public
    }
}

#[cfg(test)]
mod test {
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
        let key = LWEKey::new(params, rng);
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
        let key = crate::key::NTRUKey::new(params, rng);

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
            let key = NTRUKey::new(params, rng);
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
}

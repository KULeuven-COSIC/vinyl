//! Defining all sorts of cryptographic keys and the operations they can perform

use crate::ciphertext::{LweCiphertext, NtruVectorCiphertext};
use crate::modular::{bit_to_field, int_to_field, sample_discrete_gaussian};
use crate::params::{Params, Rng};
use crate::poly::{FFTPoly, Poly};

use swanky_field::{FiniteField, PrimeFiniteField};

/// An LWE secret key (binary)
#[derive(Debug)]
pub struct LWEKey {
    /// The actual key, stored in 64-bit chunks, from LSB to MSB
    key: Vec<u64>,

    #[cfg(debug_assertions)]
    /// If we're running in debug mode, keep track of the LWE dimension for sanity checks
    dim: usize,
}

// NOTE: this has several hardcoded assumptions about the underlying type (in the vec) being u64
// TODO: We sometimes assume the modulus fits in a single limb
impl LWEKey {
    /// Sample a new random binary LWE key
    pub fn new<BoI, BaI>(params: &Params<BoI, BaI>, rng: &mut impl Rng) -> Self {
        // We read slightly more than we really have to
        // but we can ignore the remaining parts where needed
        let mut key = vec![0; (params.dim_lwe + 63) / 64];
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
        ((self.key[i / 64] >> (i % 64)) & 1) as u8
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
        ct.b += message * int_to_field((BaI::modulus_int::<1>() >> 2).as_limbs()[0])
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

        let scale = BaI::modulus_int::<1>() >> 2;
        let half_scale = int_to_field((BaI::modulus_int::<1>() >> 3).as_limbs()[0]);

        // Add scale / 2 and floor div to round
        let out = (ct.b - mask + half_scale).into_int::<1>()
            / crypto_bigint::NonZero::new(scale).unwrap();
        out.bit_vartime(0) as u8
    }
}

/// An NTRU/NGS key, storing both the original polynomial 4f' + 1,
/// and its inverse mod X^N + 1
#[derive(Debug, Clone)]
struct NTRUKey {
    #[cfg(debug_assertions)]
    /// `4*f' + 1`, the polynomial before inversion
    f: FFTPoly,
    /// `f^-1` the actual thing used for encryption in practice
    finv: FFTPoly,
}

/// Sample a random ternary element, embedded into the modular type M
fn sample_ternary<F: PrimeFiniteField>(rng: &mut impl Rng) -> F {
    let r = rng.gen_range(0..3u8).into();
    int_to_field::<F>(r) - F::ONE
}

impl NTRUKey {
    /// Generate a fresh NTRU/NGS key of given degree
    fn new<BoI: PrimeFiniteField, BaI>(params: &Params<BoI, BaI>, rng: &mut impl Rng) -> Self {
        let mut f = Poly::<BoI>::new(1 << params.log_deg_ntru);
        let scale: BoI = int_to_field(4u8.into());

        // Until we find an invertible vector
        loop {
            f.iter_mut().for_each(|x| *x = scale * sample_ternary(rng));
            f = f + BoI::ONE;
            if let Some(finv) = f.invert() {
                return Self {
                    #[cfg(debug_assertions)]
                    f: params.fft.fwd(f),
                    finv: params.fft.fwd(finv),
                };
            }
        }
    }

    fn enc_bit_vec<BoI: PrimeFiniteField, BaI>(
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
        NtruVectorCiphertext { vec: res }
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
struct Bsk();

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
        Key {
            params,
            base,
            #[cfg(debug_assertions)]
            boot,
            public: PublicKey {
                params,
                ksk,
                bsk: Bsk(),
            },
        }
    }

    pub fn export(&self) -> &PublicKey<BoI, BaI> {
        &self.public
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::SeedableRng;

    fn rng(seed: Option<u64>) -> impl Rng {
        rand::rngs::StdRng::seed_from_u64(seed.unwrap_or(1337))
    }

    #[test]
    fn test_lwe_roundtrip() {
        let rng = &mut rng(None);
        let params = &crate::params::TESTPARAMS;
        let key = LWEKey::new(params, rng);
        for _ in 0..100 {
            let ct0 = key.encrypt(0, params, rng);
            assert_eq!(key.decrypt(ct0, params), 0);
            let ct1 = key.encrypt(1, params, rng);
            assert_eq!(key.decrypt(ct1, params), 1);
        }
    }
}

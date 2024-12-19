//! Defining all sorts of cryptographic keys and the operations they can perform

use crate::ciphertext::LweCiphertext;
use crate::modular::{bit_to_field, int_to_field, sample_discrete_gaussian};
use crate::params::{Params, Rng};
use crate::poly::Poly;

use swanky_field::{FiniteField, PrimeFiniteField, StatisticallySecureField};

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

    /// Obtain the `i`th component of the key
    /// as a bit (0 or 1 in the field).
    pub(crate) fn get<F: FiniteField>(&self, i: usize) -> F {
        #[cfg(debug_assertions)]
        debug_assert!((0..self.dim).contains(&i));

        // Entirely little endian, to the bit level
        bit_to_field(((self.key[i / 64] >> (i % 64)) & 1) as u8)
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
        // TODO: Does this impact security if we'd have a bias on the a's
        BaI: StatisticallySecureField + PrimeFiniteField,
    {
        debug_assert!(message == 0 || message == 1);

        let mut ct = LweCiphertext {
            a: Vec::with_capacity(self.dim),
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
struct NTRUKey<F> {
    f: Poly<F>,
    finv: Poly<F>,
}

/// Sample a random ternary element, embedded into the modular type M
fn sample_ternary<F: PrimeFiniteField>(rng: &mut impl Rng) -> F {
    let r = rng.gen_range(0..3u8).into();
    int_to_field::<F>(r) - F::ONE
}

impl<F: PrimeFiniteField> NTRUKey<F> {
    /// Generate a fresh NTRU/NGS key of given degree
    fn new(deg: u64, rng: &mut impl Rng) -> Self {
        let mut f = Poly::<F>::new(deg);
        let scale: F = int_to_field(4u8.into());

        // Until we find an invertible vector
        loop {
            f.iter_mut()
                .for_each(|x| *x = scale * sample_ternary(rng) + F::ONE);
            if let Some(finv) = f.invert() {
                return Self { f, finv };
            }
        }
    }
}

// TODO
// also TODO: Distinguish from LWE->LWE KSK
/// A key-switching key from NTRU to LWE
#[derive(Debug)]
struct KskNtruLwe<F>(std::marker::PhantomData<F>);

impl<F> KskNtruLwe<F> {
    /// Build a KSK that transfers from ciphertext under the ntru key
    /// to ciphertexts under the lwe key
    fn new(ntru: &NTRUKey<F>, lwe: &LWEKey) -> Self {
        todo!()
    }
}

/// A FINAL bootstrapping key: encryptions of an LWE secret key under the NTRU boot key
/// ready for use in CMUX gates.
#[derive(Debug)]
struct Bsk();

/// A full FINAL key, including all private key material
#[derive(Debug)]
pub struct Key<BootInt: 'static, BaseInt: 'static> {
    params: &'static Params<BootInt, BaseInt>,
    base: LWEKey,
    boot: NTRUKey<BootInt>,
    ksk: KskNtruLwe<BootInt>,
    bsk: Bsk,
}

impl<BoI, BaI> Key<BoI, BaI>
where
    BoI: PrimeFiniteField,
{
    pub fn new(params: &'static Params<BoI, BaI>, rng: &mut impl Rng) -> Self {
        let base = LWEKey::new(params, rng);
        let boot = NTRUKey::new(1 << params.log_deg_ntru, rng);
        let ksk = KskNtruLwe::new(&boot.clone(), &base); // Check what kind of modswitching we may need
        Key {
            params,
            base,
            boot,
            ksk,
            bsk: Bsk(),
        }
    }
}

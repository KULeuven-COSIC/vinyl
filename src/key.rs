use crate::modular::Modular;
use crate::params::Params;
use crate::poly::Poly;

// TODO: improve storage efficiency
/// An LWE secret key (binary)
#[derive(Debug)]
struct LWEKey(Vec<bool>);

impl LWEKey {
    /// Sample a new random binary LWE key
    fn new(dim: u64, rng: &mut impl random::Source) -> Self {
        let mut v = vec![false; dim as usize];
        for x in v.iter_mut() {
            *x = rng.read_u64() & 1 == 0;
        }
        Self(v)
    }
}

/// An NTRU/NGS key, storing both the original polynomial 4f' + 1,
/// and its inverse mod X^N + 1
#[derive(Debug, Clone)]
struct NTRUKey<M> {
    f: Poly<M>,
    finv: Poly<M>,
}

/// Sample a random ternary element, embedded into the modular type M
fn sample_ternary<M>(rng: &mut impl random::Source) -> M
where
    M: Modular,
    M::Base: From<u64>,
{
    loop {
        // 2^64 ≡ 1 (mod 3), so we need to exclude only the max
        let r = rng.read_u64();
        if r != u64::MAX {
            return M::construct(From::from(r % 3)) - M::construct(From::from(1));
        }
    }
}

impl<M: Modular> NTRUKey<M> {
    /// Generate a fresh NTRU/NGS key of given degree
    fn new(deg: u64, rng: &mut impl random::Source) -> Self
    where
        M::Base: From<u64> + TryInto<u64> + PartialOrd,
        M: std::fmt::Debug,
    {
        let mut f = Poly::<M>::new(deg);
        let scale = M::construct(From::from(4));
        let one = M::construct(From::from(1));

        // Until we find an invertible vector
        loop {
            f.iter_mut()
                .for_each(|x| *x = scale.clone() * sample_ternary(rng) + one.clone());
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
struct KskNtruLwe<M>(std::marker::PhantomData<M>);

impl<M: Modular> KskNtruLwe<M> {
    /// Build a KSK that transfers from ciphertext under the ntru key
    /// to ciphertexts under the lwe key
    fn new(ntru: &NTRUKey<M>, lwe: &LWEKey) -> Self {
        todo!()
    }
}

/// A FINAL bootstrapping key: encryptions of an LWE secret key under the NTRU boot key
/// ready for use in CMUX gates.
#[derive(Debug)]
struct BSK();

/// A full FINAL key, including all private key material
#[derive(Debug)]
pub struct Key<BootInt: 'static, BaseInt: 'static, NTRUExpInt: 'static> {
    // Can't use P::DIM_LWE as array size here :(
    params: &'static Params<BootInt, BaseInt, NTRUExpInt>,
    base: LWEKey,
    boot: NTRUKey<BootInt>,
    ksk: KskNtruLwe<BootInt>,
    bsk: BSK,
}

impl<BoI, BaI, NEI> Key<BoI, BaI, NEI>
where
    BoI: Modular,
    BoI: std::fmt::Debug,
    BoI::Base: From<u64> + TryInto<u64> + PartialOrd,
    BaI: Modular,
{
    pub fn new(params: &'static Params<BoI, BaI, NEI>, rng: &mut impl random::Source) -> Self {
        let base = LWEKey::new(params.dim_lwe, rng);
        let boot = NTRUKey::new(1 << params.log_deg_ntru, rng);
        let ksk = KskNtruLwe::new(&boot.clone(), &base); // Check what kind of modswitching we may need
        Key {
            params,
            base,
            boot,
            ksk,
            bsk: BSK(),
        }
    }
}

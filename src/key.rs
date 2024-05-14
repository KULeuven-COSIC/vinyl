use crate::params::{ModSwitch, Modular, Params};
use crate::poly::Poly;

#[derive(Debug)]
struct LWEKey(Vec<bool>);

impl LWEKey {
    fn new<RNG: random::Source>(dim: u64, rng: &mut RNG) -> Self {
        let mut v = vec![false; dim as usize];
        for x in v.iter_mut() {
            *x = rng.read_u64() & 1 == 0;
        }
        Self(v)
    }
}

#[derive(Debug, Clone)]
struct NTRUKey<M> {
    f: Poly<M>,
    finv: Poly<M>,
}

fn sample_ternary<M, RNG>(rng: &mut RNG) -> M
where
    M: Modular,
    M::Base: From<u64>,
    RNG: random::Source,
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
    fn new<RNG: random::Source>(deg: u64, scale: u64, rng: &mut RNG) -> Self
    where
        M::Base: From<u64> + TryInto<u64> + PartialOrd,
        M: std::fmt::Debug,
    {
        // Until we find an invertible vector
        let mut f = Poly::<M>::new(deg);
        let scale = M::construct(From::from(scale));
        let one = M::construct(From::from(1));
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
#[derive(Debug)]
struct KSK<M>(std::marker::PhantomData<M>);

impl<M: Modular> KSK<M> {
    fn new(ntru: &NTRUKey<M>, lwe: &LWEKey) -> Self {
        todo!()
    }
}

// TODO check if we can just modswitch the ntru key?
impl<M, T> ModSwitch<NTRUKey<T>> for NTRUKey<M>
where
    M: ModSwitch<T>,
{
    fn modswitch(self) -> NTRUKey<T> {
        todo!()
    }
}

#[derive(Debug)]
struct BSK();

#[derive(Debug)]
pub struct Key<BootInt: 'static, BaseInt: 'static, NTRUExpInt: 'static> {
    // Can't use P::DIM_LWE as array size here :(
    params: &'static Params<BootInt, BaseInt, NTRUExpInt>,
    base: LWEKey,
    boot: NTRUKey<BootInt>,
    ksk: KSK<BaseInt>,
    bsk: BSK,
}

impl<BoI, BaI, NEI> Key<BoI, BaI, NEI>
where
    BoI: Modular,
    BoI: std::fmt::Debug,
    BoI::Base: From<u64> + TryInto<u64> + PartialOrd,
    BaI: Modular,
{
    pub fn new<RNG: random::Source>(params: &'static Params<BoI, BaI, NEI>, rng: &mut RNG) -> Self {
        let base = LWEKey::new(params.dim_lwe, rng);
        let boot = NTRUKey::new(
            1 << params.log_deg_ntru,
            params.mod_plain * params.mod_carry,
            rng,
        );
        let ksk = KSK::new(&boot.clone().modswitch(), &base);
        Key {
            params,
            base,
            boot,
            ksk,
            bsk: BSK(),
        }
    }
}

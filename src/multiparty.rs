use swanky_field::PrimeFiniteField;

use crate::ciphertext::{LweCiphertext, MKLweCiphertext};
use crate::key::{Key, KskLweLwe, LWEKey, PublicKey};
use crate::modular::{int_to_field, ModSwitch, Modular};
use crate::params::{Params, Rng};

pub struct Client {
    key: LWEKey,
}

impl Client {
    fn new(dim: usize, rng: &mut impl Rng) -> Self {
        Self {
            key: LWEKey::new(dim, rng),
        }
    }

    pub fn encrypt<P: Params>(&self, bit: u8, rng: &mut impl Rng) -> LweCiphertext<P::BaseInt>
    where
        P::BaseInt: PrimeFiniteField,
    {
        self.key.encrypt::<P>(bit, rng)
    }
}

impl std::borrow::Borrow<LWEKey> for Client {
    fn borrow(&self) -> &LWEKey {
        &self.key
    }
}

pub struct Server<P: Params, const N: usize> {
    evaluation_key: PublicKey<P>,
    input_ksks: [KskLweLwe<P::BaseInt, 1>; N],
    output_ksk: KskLweLwe<P::BaseInt, N>,
}

impl<P: Params, const N: usize> Server<P, N>
where
    P::BaseInt: PrimeFiniteField,
    P::BaseInt: ModSwitch<P::ExpRing>,
    P::BootInt: ModSwitch<P::BaseInt>
        + PrimeFiniteField
        + for<'a> std::ops::Mul<&'a P::BootInt, Output = P::BootInt>,
    <P::ExpRing as Modular>::Single: TryInto<usize>,
    <<P::ExpRing as Modular>::Single as TryInto<usize>>::Error: std::fmt::Debug,
{
    pub fn input(&self, client: usize, ct: LweCiphertext<P::BaseInt>) -> LweCiphertext<P::BaseInt> {
        debug_assert!(client < N);

        self.input_ksks[client].keyswitch::<P>(ct)
    }

    pub fn output(&self, ct: LweCiphertext<P::BaseInt>) -> MKLweCiphertext<P::BaseInt, N> {
        self.output_ksk.keyswitch::<P>(ct)
    }

    /// bootstrap(`multiplier` * Δ/2 - (`a` + `b`))
    fn combine(
        &self,
        multiplier: impl Into<crypto_bigint::Limb>,
        a: &LweCiphertext<P::BaseInt>,
        b: &LweCiphertext<P::BaseInt>,
    ) -> LweCiphertext<P::BaseInt> {
        self.evaluation_key.bootstrap(LweCiphertext {
            a: [a.a[0]
                .iter()
                .zip(b.a[0].iter())
                .map(|(ai, bi)| *ai + *bi)
                .collect()],
            b: int_to_field::<P::BaseInt>(multiplier.into()) * P::half_scale_lwe() + a.b + b.b,
        })
    }

    pub fn not(&self, ct: &LweCiphertext<P::BaseInt>) -> LweCiphertext<P::BaseInt> {
        // Δ - ct
        LweCiphertext {
            a: [ct.a[0].iter().map(|ai| -*ai).collect()],
            b: P::scale_lwe() - ct.b,
        }
    }

    pub fn xor(
        &self,
        a: &LweCiphertext<P::BaseInt>,
        b: &LweCiphertext<P::BaseInt>,
    ) -> LweCiphertext<P::BaseInt> {
        // 2 * (a + b)
        self.evaluation_key.bootstrap(LweCiphertext {
            a: [a.a[0]
                .iter()
                .zip(b.a[0].iter())
                .map(|(ai, bi)| *ai + *bi + *ai + *bi)
                .collect()],
            b: a.b + b.b + a.b + b.b,
        })
    }

    pub fn nand(
        &self,
        a: &LweCiphertext<P::BaseInt>,
        b: &LweCiphertext<P::BaseInt>,
    ) -> LweCiphertext<P::BaseInt> {
        self.combine(3u64, a, b)
    }

    pub fn and(
        &self,
        a: &LweCiphertext<P::BaseInt>,
        b: &LweCiphertext<P::BaseInt>,
    ) -> LweCiphertext<P::BaseInt> {
        self.combine(7u64, a, b)
    }

    pub fn or(
        &self,
        a: &LweCiphertext<P::BaseInt>,
        b: &LweCiphertext<P::BaseInt>,
    ) -> LweCiphertext<P::BaseInt> {
        self.combine(1u64, a, b)
    }
}

pub fn setup<P: Params, const N: usize>(rng: &mut impl Rng) -> ([Client; N], Server<P, N>)
where
    P::BootInt: PrimeFiniteField,
    P::BaseInt: PrimeFiniteField,
{
    let clients = std::array::from_fn(|_| Client::new(P::DIM_LWE, rng));
    let fhe_key = Key::new(rng);
    // TODO: improve on these clones
    let input_ksks =
        std::array::from_fn(|i| KskLweLwe::new::<P, _>(&clients[i].key, &[&fhe_key.base], rng));
    let output_ksk = KskLweLwe::new::<P, _>(
        &fhe_key.base,
        &std::array::from_fn(|i| &clients[i].key),
        rng,
    );
    let server = Server {
        evaluation_key: fhe_key.take_public(),
        input_ksks,
        output_ksk,
    };

    (clients, server)
}

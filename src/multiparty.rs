use swanky_field::PrimeFiniteField;

use crate::ciphertext::{LweCiphertext, MKLweCiphertext};
use crate::key::{Key, KskLweLwe, KskNoise, KskNtruMKLwe, LWEKey, PublicKey};
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
    // output_ksk: KskLweLwe<P::BaseInt, N>,
    output_ksk: KskNtruMKLwe<P::BootInt, N>,
    secret_key: LWEKey, // TODO
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
        // self.output_ksk.keyswitch::<P>(ct)
        let (a, b) = ct.unpack();
        self.evaluation_key.bootstrap_with_ksk(
            LweCiphertext {
                a: [a.into_iter().map(|ai| ai + ai).collect()],
                b: b + b,
            },
            &self.output_ksk,
        )
    }

    pub fn output_noise(
        &self,
        ct: LweCiphertext<P::BaseInt>,
        m: u8,
        clients: &[Client; N],
    ) -> (usize, usize) {
        fn noise<P: Params, const N: usize>(
            ct: MKLweCiphertext<P::BaseInt, N>,
            key: &[impl std::borrow::Borrow<LWEKey>; N],
            m: u8,
        ) -> usize
        where
            P::BaseInt: PrimeFiniteField,
        {
            let mask =
                ct.a.iter()
                    .zip(key.iter())
                    .flat_map(|(as_, k)| {
                        as_.iter()
                            .zip(k.borrow().iter(P::DIM_LWE))
                            .map(|(ai, si): (&P::BaseInt, P::BaseInt)| *ai * si)
                    })
                    .sum::<P::BaseInt>();
            let err = crate::modular::lift_centered(
                ct.b + mask - P::scale_lwe() * int_to_field(m.into()),
            );
            if err.abs() == 0 {
                0
            } else {
                err.abs().ilog2() as usize
            }
        }

        (
            noise::<P, 1>(ct.clone(), &[&self.secret_key], m),
            noise::<P, N>(self.output(ct), clients, m),
        )
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
    let mut clients = std::array::from_fn(|_| Client::new(P::DIM_LWE, rng));
    let (fhe_key, ntru_key) = Key::new_and_ntru(rng);
    clients[0] = Client {
        key: fhe_key.base.clone(),
    };
    // TODO: improve on these clones
    let input_ksks = std::array::from_fn(|i| {
        KskLweLwe::new::<P, _>(&clients[i].key, &[&fhe_key.base], rng, KskNoise::PerParty)
    });
    // let output_ksk = KskLweLwe::new::<P, _>(
    // &fhe_key.base,
    // &std::array::from_fn(|i| &clients[i].key),
    // rng,
    // KskNoise::Single,
    // );
    let output_ksk = KskNtruMKLwe::new_mk::<P>(&ntru_key, &clients, rng, KskNoise::Single);
    let secret_key = fhe_key.base.clone(); // TODO
    let server = Server {
        evaluation_key: fhe_key.take_public(),
        input_ksks,
        output_ksk,
        secret_key, // TODO
    };

    (clients, server)
}

use std::time::Duration;

use rand::{rngs::OsRng, Rng};
use vinyl::{
    bench::{fft::FFTPlan, poly::Poly},
    ciphertext::NtruScalarCiphertext,
    key::{Bsk, Key, KskLweLwe, KskNoise, KskNtruLwe, LWEKey, NTRUKey},
    params::{FinalParams, Params},
};
type BootInt = <FinalParams as Params>::BootInt;
use swanky_field::FiniteRing;

use criterion::{criterion_group, criterion_main, Criterion};

criterion_main!(
    bench_lwe,
    bench_fft,
    bench_bootstrap,
    bench_keygen,
    bench_multiparty
);

criterion_group!(bench_lwe, lwe_enc, lwe_dec);
fn lwe_enc(c: &mut Criterion) {
    let rng = &mut OsRng;
    let setup_rng = &mut OsRng;
    let key = LWEKey::new(FinalParams::DIM_LWE, rng);

    c.bench_function("lwe_enc", |bencher| {
        bencher.iter_batched(
            || setup_rng.gen_range(0..=1),
            |b| key.encrypt::<FinalParams>(b, rng),
            criterion::BatchSize::SmallInput,
        );
    });
}

fn lwe_dec(c: &mut Criterion) {
    let rng = &mut OsRng;
    let key = LWEKey::new(FinalParams::DIM_LWE, rng);

    c.bench_function("lwe_dec", |bencher| {
        bencher.iter_batched(
            || key.encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
            |ct| key.decrypt::<FinalParams>(ct),
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(bench_fft, fft_fwd, fft_inv, complex_mult);
fn fft_fwd(c: &mut Criterion) {
    let plan = FFTPlan::new(1 << FinalParams::LOG_DEG_NTRU);
    let rng = &mut OsRng;

    c.bench_function("fft_fwd", |bencher| {
        bencher.iter_batched(
            || Poly::<BootInt>::rand(1 << FinalParams::LOG_DEG_NTRU, rng),
            |poly| plan.fwd(poly),
            criterion::BatchSize::SmallInput,
        )
    });
}

fn fft_inv(c: &mut Criterion) {
    let plan = FFTPlan::new(1 << FinalParams::LOG_DEG_NTRU);
    let rng = &mut OsRng;

    c.bench_function("fft_inv", |bencher| {
        bencher.iter_batched(
            || plan.fwd(Poly::<BootInt>::rand(1 << FinalParams::LOG_DEG_NTRU, rng)),
            |poly| plan.inv::<BootInt>(poly),
            criterion::BatchSize::SmallInput,
        )
    });
}

fn complex_mult(c: &mut Criterion) {
    let plan = FFTPlan::new(1 << FinalParams::LOG_DEG_NTRU);
    let rng = &mut OsRng;

    c.bench_function("complex_mult", |bencher| {
        bencher.iter_batched(
            || {
                (
                    plan.fwd(Poly::<BootInt>::rand(1 << FinalParams::LOG_DEG_NTRU, rng)),
                    plan.fwd(Poly::<BootInt>::rand(1 << FinalParams::LOG_DEG_NTRU, rng)),
                )
            },
            |(a, b)| a * b,
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(
    name = bench_bootstrap;
    config = Criterion::default().measurement_time(Duration::from_secs(20));
    targets = single_bootstrap, double_bootstrap, /* ext_prod, */ keyswitch, keyswitch_lwe
);

fn single_bootstrap(c: &mut Criterion) {
    let rng = &mut OsRng;
    let private = Key::<FinalParams>::new(rng);
    let public = private.export();

    c.bench_function("single_bootstrap", |bencher| {
        bencher.iter_batched(
            || private.encrypt(rng.gen_range(0..=1), rng),
            |ct| public.bootstrap(ct),
            criterion::BatchSize::SmallInput,
        )
    });
}

fn double_bootstrap(c: &mut Criterion) {
    let rng = &mut OsRng;
    let private = Key::<FinalParams>::new(rng);
    let public = private.export();

    c.bench_function("double_bootstrap", |bencher| {
        bencher.iter_batched(
            || private.encrypt(rng.gen_range(0..=1), rng),
            |ct| public.bootstrap(public.bootstrap(ct)),
            criterion::BatchSize::SmallInput,
        )
    });
}

// fn ext_prod(c: &mut Criterion) {
//     let rng = &mut OsRng;
//     let len = 1 << FinalParams::LOG_DEG_NTRU;
//     let (key, _) = NTRUKey::new::<FinalParams>(rng);
//     let fft = FinalParams::fft();

//     c.bench_function("external_product", |bencher| {
//         bencher.iter_batched(
//             || {
//                 let mut scalar = Poly::<BootInt>::new(len);
//                 scalar.iter_mut().for_each(|x| *x = sample_ternary(rng));
//                 let x = NtruScalarCiphertext::trivial::<FinalParams>(scalar.clone());
//                 let mut vector = Poly::<BootInt>::new(len);
//                 vector.0[rng.gen_range(0..len)] = sample_ternary(rng);
//                 let y = NtruVectorCiphertext::trivial::<FinalParams>(vector.clone());
//             },
//             |()| (),
//             criterion::BatchSize::SmallInput,
//         )
//     });
// }

fn keyswitch(c: &mut Criterion) {
    let rng = &mut OsRng;
    let lwe = LWEKey::new(FinalParams::DIM_LWE, rng);
    let (ntru, ntru_coefs) = NTRUKey::new::<FinalParams>(rng);
    let scalar_one = NtruScalarCiphertext::trivial::<FinalParams>(
        Poly::new(1 << FinalParams::LOG_DEG_NTRU) + BootInt::ONE,
    );
    let ksk = KskNtruLwe::new::<FinalParams>(&ntru_coefs, &lwe, rng);
    let fft = FinalParams::fft();

    c.bench_function("keyswitch (ngs->lwe)", |bencher| {
        bencher.iter_batched(
            || {
                scalar_one
                    .clone()
                    .external_product::<FinalParams>(&ntru.enc_bit_vec::<FinalParams>(0, rng), fft)
            },
            |ct| ksk.key_switch::<FinalParams>(ct),
            criterion::BatchSize::SmallInput,
        )
    });
}

fn keyswitch_lwe(c: &mut Criterion) {
    let rng = &mut OsRng;
    let lwe1 = LWEKey::new(FinalParams::DIM_LWE, rng);
    let lwe2 = LWEKey::new(FinalParams::DIM_LWE, rng);
    let ksk = KskLweLwe::new::<FinalParams, _>(&lwe1, &[lwe2], rng, KskNoise::Single);

    c.bench_function("keyswitch (lwe->lwe)", |bencher| {
        bencher.iter_batched(
            || lwe1.encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
            |ct| ksk.keyswitch::<FinalParams>(ct),
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(bench_keygen, lwegen, ntrugen, bskgen, ntrulwekskgen, keygen);

fn lwegen(c: &mut Criterion) {
    let rng = &mut OsRng;
    c.bench_function("lwegen", |bencher| {
        bencher.iter(|| LWEKey::new(FinalParams::DIM_LWE, rng))
    });
}

fn ntrugen(c: &mut Criterion) {
    let rng = &mut OsRng;
    c.bench_function("ntrugen", |bencher| {
        bencher.iter(|| NTRUKey::new::<FinalParams>(rng))
    });
}

fn bskgen(c: &mut Criterion) {
    let rng = &mut OsRng;
    let lwe = LWEKey::new(FinalParams::DIM_LWE, rng);
    let (ntru, _) = NTRUKey::new::<FinalParams>(rng);

    c.bench_function("bskgen", |bencher| {
        bencher.iter(|| Bsk::new::<FinalParams>(&ntru, &lwe, rng))
    });
}

fn ntrulwekskgen(c: &mut Criterion) {
    let rng = &mut OsRng;
    let lwe = LWEKey::new(FinalParams::DIM_LWE, rng);
    let (_, ntru_coefs) = NTRUKey::new::<FinalParams>(rng);

    c.bench_function("ntrulwekskgen", |bencher| {
        bencher.iter(|| KskNtruLwe::new::<FinalParams>(&ntru_coefs, &lwe, rng))
    });
}

fn keygen(c: &mut Criterion) {
    let rng = &mut OsRng;

    c.bench_function("keygen", |bencher| {
        bencher.iter(|| Key::<FinalParams>::new(rng))
    });
}

criterion_group!(
    bench_multiparty,
    multiparty_keygen,
    multiparty_input,
    multiparty_output,
    multiparty_mk_decrypt,
    multiparty_gates
);

macro_rules! multi {
    ($b:ident) => {
        $b!(2);
        $b!(4);
        $b!(8);
        $b!(16);
        $b!(32);
    };
}

fn multiparty_keygen(c: &mut Criterion) {
    macro_rules! b {
        ($n:expr) => {
            let rng = &mut OsRng;
            c.bench_function(&format!("multiparty_keygen ({} parties)", $n), |bencher| {
                bencher.iter(|| vinyl::multiparty::setup::<FinalParams, $n>(rng))
            });
        };
    }
    multi!(b);
}

fn multiparty_input(c: &mut Criterion) {
    macro_rules! b {
        ($n:expr) => {
            let rng = &mut OsRng;
            let (clients, server) = vinyl::multiparty::setup::<FinalParams, $n>(rng);
            c.bench_function(
                &format!("multiparty input ({} parties)", $n),
                move |bencher| {
                    bencher.iter_batched(
                        || {
                            let i = rng.gen_range(0..$n);
                            (
                                i,
                                clients[i].encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
                            )
                        },
                        |(i, ct)| {
                            server.input(i, ct);
                        },
                        criterion::BatchSize::LargeInput,
                    )
                },
            );
        };
    }
    multi!(b);
}

fn multiparty_output(c: &mut Criterion) {
    macro_rules! b {
        ($n:expr) => {
            let rng = &mut OsRng;
            let (clients, server) = vinyl::multiparty::setup::<FinalParams, $n>(rng);
            c.bench_function(
                &format!("multiparty output ({} parties)", $n),
                move |bencher| {
                    bencher.iter_batched(
                        || {
                            let i = rng.gen_range(0..$n);
                            server.input(
                                i,
                                clients[i].encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
                            )
                        },
                        |ct| {
                            server.output(ct);
                        },
                        criterion::BatchSize::LargeInput,
                    )
                },
            );
        };
    }
    multi!(b);
}

fn multiparty_mk_decrypt(c: &mut Criterion) {
    macro_rules! b {
        ($n:expr) => {
            let rng = &mut OsRng;
            let (clients, server) = vinyl::multiparty::setup::<FinalParams, $n>(rng);
            c.bench_function(
                &format!("multiparty decrypt ({} parties)", $n),
                move |bencher| {
                    bencher.iter_batched(
                        || {
                            let i = rng.gen_range(0..$n);
                            server.output(server.input(
                                i,
                                clients[i].encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
                            ))
                        },
                        |ct| {
                            ct.decrypt::<FinalParams, _>(&clients);
                        },
                        criterion::BatchSize::SmallInput,
                    )
                },
            );
        };
    }
    multi!(b);
}

fn multiparty_gates(c: &mut Criterion) {
    macro_rules! b {
        ($n:expr) => {
            let rng = &mut OsRng;
            let (clients, server) = vinyl::multiparty::setup::<FinalParams, $n>(rng);
            c.bench_function(
                &format!("multiparty nand gate ({} parties)", $n),
                move |bencher| {
                    bencher.iter_batched(
                        || {
                            let i = rng.gen_range(0..$n);
                            let j = rng.gen_range(0..$n);
                            (
                                server.input(
                                    i,
                                    clients[i].encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
                                ),
                                server.input(
                                    j,
                                    clients[j].encrypt::<FinalParams>(rng.gen_range(0..=1), rng),
                                ),
                            )
                        },
                        |(a, b)| server.nand(&a, &b),
                        criterion::BatchSize::SmallInput,
                    )
                },
            );
        };
    }
    multi!(b);
}

// TODO? Look at performance of some of the modular stuff

use rand::{rngs::OsRng, Rng};
use vinyl::{
    bench::{fft::FFTPlan, poly::Poly},
    key::{Bsk, Key, KskNtruLwe, LWEKey, NTRUKey},
    params::{FinalParams, Params},
};
type BootInt = <FinalParams as Params>::BootInt;

use criterion::{criterion_group, criterion_main, Criterion};

criterion_main!(bench_lwe, bench_fft, bench_bootstrap);

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

criterion_group!(bench_bootstrap, single_bootstrap, double_bootstrap);

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

// TODO? Look at performance of some of the modular stuff

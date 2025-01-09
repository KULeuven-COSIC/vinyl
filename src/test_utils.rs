#![cfg(test)]

use std::ops::{AddAssign, Mul, Sub};

use crate::params::Rng;

use rand::SeedableRng;

pub(crate) fn rng(seed: Option<u64>) -> impl Rng {
    rand::rngs::StdRng::seed_from_u64(seed.unwrap_or(1337))
}

pub(crate) fn slow_negacyclic_mult<
    F: Clone + Copy + AddAssign<F> + Mul<F, Output = F> + Sub<F, Output = F>,
>(
    a: &[F],
    b: &[F],
    zero: F,
) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    let n = a.len();

    let mut full = vec![zero; 2 * n];
    for i in 0..n {
        for j in 0..n {
            full[i + j] += a[i] * b[j];
        }
    }

    let mut res = Vec::with_capacity(n);
    for i in 0..n {
        res.push(full[i] - full[i + n]);
    }
    res
}

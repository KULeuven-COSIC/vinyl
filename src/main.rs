use rand::SeedableRng;
use vinyl::key::Key;
use vinyl::params::*;

fn main() {
    let key = Key::new(&TESTPARAMS, &mut rand::rngs::StdRng::seed_from_u64(42));
}

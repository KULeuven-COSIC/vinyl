use vinyl::key::Key;
use vinyl::params::*;

fn main() {
    let key = Key::new(&TESTPARAMS, &mut random::default(42));
}

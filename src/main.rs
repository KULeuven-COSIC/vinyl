use vinyl::key::Key;
use vinyl::params::*;

fn main() {
    let key = Key::new(&MESSAGE_7_CARRY_0, &mut random::default(42));
}

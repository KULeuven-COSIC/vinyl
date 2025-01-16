use vinyl::key::Key;
use vinyl::params::*;

fn main() {
    let rng = &mut rand::rngs::OsRng;
    let private = Key::<FinalParams>::new(rng);
    let public = private.export();

    let ct = private.encrypt(0, rng);
    let bootstrapped_ct = public.bootstrap(ct);
    let bootstrapped_plain = private.decrypt(bootstrapped_ct);

    println!(
        "Encrypted, bootstrapped and decrypted back to {}",
        bootstrapped_plain
    );
}

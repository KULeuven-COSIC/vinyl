use vinyl::params::*;

fn main() {
    let rng = &mut rand::rngs::OsRng;
    type Params = FinalParams;
    let (clients, server) = vinyl::multiparty::setup::<Params, 2>(rng);
    for b0 in 0..=1 {
        let ct0 = server.input(0, clients[0].encrypt::<Params>(b0, rng));
        for b1 in 0..=1 {
            let ct1 = server.input(1, clients[1].encrypt::<Params>(b1, rng));

            let not0 = server.not(&ct0);
            let not1 = server.not(&ct1);
            let xor = server.xor(&ct0, &ct1);
            let nand = server.nand(&ct0, &ct1);
            let and = server.and(&ct0, &ct1);
            let or = server.or(&ct0, &ct1);

            for (ct, v) in [
                (not0, 1 - b0),
                (not1, 1 - b1),
                (xor, b0 ^ b1),
                (nand, 1 - (b0 & b1)),
                (and, b0 & b1),
                (or, b0 | b1),
            ] {
                eprintln!(".");
                assert_eq!(server.output(ct).decrypt::<Params, _>(&clients), v);
            }
        }
    }
}

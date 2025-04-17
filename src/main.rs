use vinyl::params::*;

fn main() {
    let rng = &mut rand::rngs::OsRng;
    type Params = FinalParams;
    let (clients, server) = vinyl::multiparty::setup::<Params, 2>(rng);
    println!("Ready");
    for line in std::io::stdin().lines() {
        let line = line.unwrap();
        let spl = line.split_whitespace().collect::<Vec<_>>();
        let b0 = server.input(
            0,
            clients[0].encrypt::<Params>(spl[1].parse().unwrap(), rng),
        );
        let ct = if spl[0] == "not" {
            server.not(&b0)
        } else {
            let b1 = server.input(
                1,
                clients[1].encrypt::<Params>(spl[2].parse().unwrap(), rng),
            );
            if spl[0] == "xor" {
                server.xor(&b0, &b1)
            } else if spl[0] == "nand" {
                server.nand(&b0, &b1)
            } else if spl[0] == "and" {
                server.and(&b0, &b1)
            } else if spl[0] == "or" {
                server.or(&b0, &b1)
            } else {
                println!("Unknown gate {}", spl[0]);
                continue;
            }
        };
        println!("~> {}", server.output(ct).decrypt::<Params, _>(&clients));
    }
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
                assert_eq!(server.output(ct).decrypt::<Params, _>(&clients), v);
            }
        }
    }
}

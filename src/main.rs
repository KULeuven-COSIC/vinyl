use rand::Rng;
use vinyl::params::*;

struct Measure(String, std::time::Instant);
impl Measure {
    fn start(label: impl Into<String>) -> Self {
        Self(label.into(), std::time::Instant::now())
    }
}
impl Drop for Measure {
    fn drop(&mut self) {
        println!("{} took {:.2?}", self.0, self.1.elapsed());
    }
}

macro_rules! measure {
    ($label:expr, $rest:expr) => {{
        let _m = Measure::start($label);
        $rest
    }};
}

fn main() {
    let rng = &mut rand::rngs::OsRng;
    type Params = FinalParams;
    const N: usize = 8;
    println!("Number of parties: {N}");
    let (clients, server) = measure!(
        "Multiparty setup",
        vinyl::multiparty::setup::<Params, N>(rng)
    );

    let mut bits = Vec::with_capacity(N * 128);
    {
        let _m = Measure::start("Input 128 bits per party");
        for i in 0..128 {
            for p in 0..N {
                bits.push(server.input(p, clients[p].encrypt::<Params>(i & 1, rng)));
            }
        }
    }

    {
        let _m = Measure::start("100000 random gates");
        for _ in 0..100000 {
            let a = rng.gen_range(0..(N * 128));
            let b = rng.gen_range(0..(N * 128));
            bits[a] = server.nand(&bits[a], &bits[b]);
        }
    }

    let mut outputs = Vec::with_capacity(N * 128);
    {
        let _m = Measure::start("Output bits");
        for b in bits {
            outputs.push(server.output(b));
        }
    }

    let _ = measure!(
        "Decrypt outputs",
        outputs.into_iter().for_each(|o| {
            o.decrypt::<Params, _>(&clients);
        })
    );
}

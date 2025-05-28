# Vinyl <sup><sub>(*/vī′nəl/*)</sub></sup>

A toy implementation of single-key and multiparty FHE, based on [FINAL](https://doi.org/10.1007/978-3-031-22966-4_7),
accompanying the paper [Multi-Party FHE Redefined: A Framework for Unlimited Participants](https://eprint.iacr.org/2025/965).

## Running

Ensure a rust compiler is available, along with cargo (tested with version `1.85.0`).
The `m4` preprocessor and `make` are also necessary for the dependencies to build.
- To run an example program that performs some multiparty FHE operations and outputs some timing information: `cargo run --release`
- To run the tests: `cargo test --release` (without the `--release` flag, the tests will also work and even check a few extra things, but take significantly longer).
- To run the benchmarks: `cargo bench -Fbench` (the `bench` feature changes the visibility of some functions so that some more granular benchmarks can inspect the performance of smaller components). Afterwards, you can additionally open `./target/criterion/report/index.html` in a browser to see the benchmark report as a web page.

mod variant_calling;

use std::env;
use variant_calling::call_variants;

fn main() {
    let args: Vec<String> = env::args().collect();
    let args: Vec<&str> = args.iter().map(|x| &x[..]).collect();
    let bams = &args[1..];
    call_variants(&bams);
}

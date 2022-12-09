/// Quantum Chemistry program.
///
/// Authors:
/// - Dmytruk Mykhailo
/// - Kuznetzov Volodymyr
/// - Klyuyev  Illia
/// TODO: write it in `Cargo.toml` file.

#[macro_use]

// use std::{
//     env,
//     fs::read_to_string as read_file,
// };

mod atom;
mod calculations;
mod eigh;
mod math;
mod molecule;
mod utils;
mod vec3d;
use crate::{
    // atom::Atom,
    calculations::{run_hf, compare},
    math::float,
    molecule::Molecule,
    // utils::ExtensionUnwrapOrExit,
    // vec3d::Vec3d,
};


fn main() {
    println!("
\n
Welcome to QuCy - Quantum Chemistry program.

Authors:
- Dmytruk Mykhailo
- Kuznetzov Volodymyr
- Klyuyev Illia
\n
    ");
    // let args: Vec<String> = env::args().collect();
    // let file_path = args.get(1).unwrap_or_exit_with_msg("No file name given");
    // let contents = read_file(file_path).unwrap_or_exit_with_msg("Unable to open specified file");
    // // println!("File contents:\n{contents}");
    // let molecule: Molecule = Molecule::from_string(&contents);
    // dbg!(molecule);


    // Test H2:
    test("
0           # total charge
H 0 0 0     # first  atom and it's xyz
H 0 0 1.4   # second atom and it's xyz
        ",
        -1.11675930740
    );

    // Test HeH+
    test("
1
He 0 0 0
H  0 0 1.4632
        ",
        -2.8606621637
    );
    // Test HeH+
    test("
0
O -0.5 0 0
O +0.5 0 0
        ",
        -149.4617
    );
}

fn test(input_text: &str, hf_e_ref: float) {
    println!("Molecule:{}", input_text);
    let mol = Molecule::from_string(input_text);
    let hf_e = run_hf(&mol);
    compare(hf_e, hf_e_ref, None);
}


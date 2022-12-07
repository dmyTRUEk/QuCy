/// Quantum Chemistry program.
///
/// Authors:
/// - Dmytruk Mykhailo
/// - Kuznetzov Volodymyr
/// - Kluev Illia
/// TODO: write it in `Cargo.toml` file.

// use std::{
//     env,
//     fs::read_to_string as read_file,
// };

mod atom;
mod calculations;
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
    // let args: Vec<String> = env::args().collect();
    // let file_path = args.get(1).unwrap_or_exit_with_msg("No file name given");
    // let contents = read_file(file_path).unwrap_or_exit_with_msg("Unable to open specified file");
    // // println!("File contents:\n{contents}");
    // let molecule: Molecule = Molecule::from_string(&contents);
    // dbg!(molecule);
    test1();
}

/// Test of H2
fn test1() {
    let input_text: &str = "
        0
        H 0 0 0
        H 0 0 1.4
    ";
    let mol = Molecule::from_string(input_text);
    let hf_e = run_hf(&mol);
    const hf_e_ref: float = -1.11675930740;
    compare(hf_e, hf_e_ref, None);
}


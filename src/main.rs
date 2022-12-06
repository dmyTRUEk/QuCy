/// Quantum Chemistry program.
///
/// Authors:
/// - Dmytruk Mykhailo
/// - Kuznetzov Volodymyr
/// - Kluev Illia
/// TODO: write it in `Cargo.toml` file.

use std::{
    env,
    fs::read_to_string as read_file,
};

mod atom;
mod math;
mod utils;
mod vec3d;
use crate::{
    atom::Atom,
    utils::ExtensionUnwrapOrExit,
    vec3d::Vec3d,
};


fn main() {
    let args: Vec<String> = env::args().collect();
    let file_path = args.get(1).unwrap_or_exit_with_msg("No file name given");
    let contents = read_file(file_path).unwrap_or_exit_with_msg("Unable to open specified file");
    // println!("File contents:\n{contents}");
    let mut atoms: Vec<Atom> = vec![];
    for line in contents.lines() {
        let content: String = line.split('#').collect::<Vec<&str>>().get(0).unwrap().trim().to_string();
        if content.is_empty() { continue; }
        // println!("{content}");
        let parts: Vec<&str> = content.split(' ').collect();
        let atom: Atom = match parts.len() {
            0 => { continue; }
            1 => { // name
                if !atoms.is_empty() { panic!("Only first atom without coorinates allowed.") }
                Atom::new(
                    parts[0].to_string(),
                    Vec3d::new(0.0, 0.0, 0.0),
                    None
                )
            }
            4 => { // name, x, y, z
                Atom::new(
                    parts[0].to_string(),
                    Vec3d::new(
                        parts[1].parse::<f64>().unwrap(),
                        parts[2].parse::<f64>().unwrap(),
                        parts[3].parse::<f64>().unwrap(),
                    ),
                    None
                )
            }
            5 => { // name, mass, x, y, z
                todo!()
            }
            _ => { panic!("parts.len() = {l}, parts = {parts:?}", l=parts.len()) }
        };
        atoms.push(atom);
    }

    println!("{atoms:#?}");
}


/// Molecule

use crate::{
    atom::Atom,
    calculations::CGF,
    math::float,
    vec3d::Vec3d,
};


#[derive(Debug, Clone)]
pub struct Molecule {
    // TODO: refactor to getters
    pub atoms: Vec<Atom>,
    // number of atoms == `atoms.len()`
    pub cgfs: Vec<CGF>,
    // number_of_electrons: usize,
}
impl Molecule {
    pub fn from_string(string: &str) -> Self {
        let (atoms, self_charge): (Vec<Atom>, i32) =
            Molecule::parse_string_into_data(string);
        let total_atoms_charge: float = atoms.iter().map(|a| a.get_charge()).sum();
        assert!(total_atoms_charge.fract() < 1e-3);
        let total_atoms_charge: i32 = total_atoms_charge.round() as i32;
        // let number_of_electrons: i32 = total_atoms_charge - self_charge;
        // assert!(0 <= number_of_electrons);
        // assert!((number_of_electrons as u128) < usize::MAX as u128);
        // TODO: remove asserts and use `.into()` or `usize::from(...)`
        // let number_of_electrons: usize = number_of_electrons as usize;
        let mut cgfs: Vec<CGF> = vec![];
        for atom in &atoms {
            for z in atom.get_zetas() {
                cgfs.push(CGF::from(
                    z,
                    atom.get_priciple_quantum_number(),
                    atom.get_position(),
                ));
            }
        }
        Self {
            atoms,
            cgfs,
            // number_of_electrons,
        }
    }

    fn parse_string_into_data(string: &str) -> (Vec<Atom>, i32) {
        let mut atoms: Vec<Atom> = vec![];
        let mut self_charge: Option<i32> = None;
        for line in string.lines() {
            let content: String = line.split('#').collect::<Vec<&str>>().get(0).unwrap().trim().to_string();
            if content.is_empty() { continue; }
            // println!("{content}");
            let parts: Vec<&str> = content.split(' ').collect();
            let atom: Atom = match parts.len() {
                0 => { continue; }
                1 => { // just name or self_charge
                    if !atoms.is_empty() { panic!("Only first atom without coorinates allowed.") }
                    if let Ok(self_charge_) = parts[0].to_string().parse::<i32>() {
                        self_charge = Some(self_charge_);
                        continue;
                    }
                    else {
                        Atom::new(
                            parts[0].to_string(),
                            Vec3d::new(0.0, 0.0, 0.0),
                            None
                        )
                    }
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
        // println!("{atoms:#?}");
        const DEFAULT_MOLECULE_CHARGE: i32 = 0; // = - number of additional electrons
        (
            atoms,
            self_charge.unwrap_or(DEFAULT_MOLECULE_CHARGE)
        )
    }

    pub fn get_number_of_atoms(&self) -> usize { self.atoms.len() }
}


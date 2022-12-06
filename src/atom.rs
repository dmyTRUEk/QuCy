/// Atom

use crate::vec3d::Vec3d;


// TODO?: AtomType enum


#[derive(Debug, Clone)]
pub struct Atom {
    name: String,
    mass: f64,
    position: Vec3d,
}
impl Atom {
    pub fn new(
        name: String,
        position: Vec3d,
        mass: Option<f64>,
    ) -> Self {
        Self {
            position,
            name: name.clone(),
            mass: mass.unwrap_or(Atom::name_to_mass(&name)),
        }
    }

    pub fn name_to_mass(name: &str) -> f64 {
        match name {
            "H"  => { 1.0 }
            "He" => { 2.0 }
            "O"  => { 8.0 }
            _ => { todo!() }
        }
    }
}


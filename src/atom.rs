/// Atom related stuff

use crate::{
    math::float,
    vec3d::Vec3d,
};


#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomType {
    H,                                                                 He,
    Li, Be,                                        B , C , N , O , F , Ne,
    Na, Mg,                                        Al, Si, P , S , Cl, Ar,
    K , Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
}


#[derive(Debug, Clone, Copy)]
pub struct Atom {
    atom_type: AtomType,
    // mass: float,
    position: Vec3d,
    // zeta: float,
    // charge: float,
    // principle_quantum_number: float,
}
impl Atom {
    pub fn new(
        name: String,
        position: Vec3d,
        // mass: Option<float>,
    ) -> Self {
        let atom_type: AtomType = Atom::name_to_type(&name);
        // let mass: float = mass.unwrap_or(Atom::atom_type_to_mass(atom_type));
        Self { atom_type, position }
    }

    // TODO?: make special type for `AtomMass=float`, etc
    // pub const fn get_type(&self) -> AtomType { self.atom_type }
    // pub const fn get_mass(&self) -> float { self.mass }
    pub const fn get_position(&self) -> Vec3d { self.position }
    pub fn get_zetas(&self) -> Vec<float> { Self::atom_type_to_zeta(self.atom_type) }
    pub const fn get_charge(&self) -> float { Self::atom_type_to_charge(self.atom_type) }
    pub const fn get_priciple_quantum_number(&self) -> i32 { Self::atom_type_to_principle_quantum_number(self.atom_type) }

    pub fn name_to_type(name: &str) -> AtomType {
        match name {
            "H"  => { AtomType::H }
            "He" => { AtomType::He }
            "Li" => { AtomType::Li }
            "O"  => { AtomType::O }
            _ => { todo!() }
        }
    }

    // pub const fn atom_type_to_mass(atom_type: AtomType) -> float {
    //     match atom_type {
    //         AtomType::H  => { 1.0 }
    //         AtomType::He => { 2.0 }
    //         AtomType::Li => { 3.0 }
    //         AtomType::O  => { 8.0 }
    //         _ => { todo!() }
    //     }
    // }

    pub fn atom_type_to_zeta(atom_type: AtomType) -> Vec<float> {
        match atom_type {
            AtomType::H  => { vec![1.24] }
            AtomType::He => { vec![2.0925] }
            AtomType::Li => { todo!() }
            AtomType::O  => { todo!() }
            _ => { todo!() }
        }
    }

    pub const fn atom_type_to_charge(atom_type: AtomType) -> float {
        match atom_type {
            AtomType::H  => { 1.0 }
            AtomType::He => { 2.0 }
            AtomType::Li => { 3.0 }
            AtomType::O  => { 8.0 }
            _ => { todo!() }
        }
    }

    pub const fn atom_type_to_principle_quantum_number(atom_type: AtomType) -> i32 {
        match atom_type {
            AtomType::H  => { 1 }
            AtomType::He => { 1 }
            AtomType::Li => { todo!() }
            AtomType::O  => { todo!() }
            _ => { todo!() }
        }
    }
}


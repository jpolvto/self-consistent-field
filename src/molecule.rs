use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    pub atomic_number: i32,
    pub position: Vector3<f32>,
    number_of_electrons: i32,
    principal_quantum_number: i32
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>
}

impl Molecule {
    fn get_nuclear_repulsion_energy(&self) -> f32 {
        let mut energy: f32 = 0.0;

        for atom1 in 0..self.atoms.len() {
            for atom2 in (atom1 + 1)..self.atoms.len() {
                energy += (self.atoms[atom1].atomic_number) as f32 / (self.atoms[atom1].position - self.atoms[atom2].position).norm();
            } 
        }
        return energy;
    }
}

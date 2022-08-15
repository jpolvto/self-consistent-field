use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::{basisset::BasisSet, gaussian::{SlaterOrbital}};

#[derive(Serialize, Deserialize, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    // every atom must have an atomic number to identify the atom
    pub atomic_number: i32,

    //the position is specified as a 3D vector
    pub position: Vector3<f32>,

    //open shell calculations
    number_of_electrons: i32,

    //open shell calculations
    principal_quantum_number: i32
}

impl Molecule {
    fn get_nuclear_repulsion_energy(&self) -> f32 {
        let mut energy: f32 = 0.0;

        for atom_a in 0..self.atoms.len() {
            for atom_b in (atom_a + 1)..self.atoms.len() {
                energy += (self.atoms[atom_a].atomic_number) as f32 / (self.atoms[atom_a].position - self.atoms[atom_b].position).norm();
            } 
        }
        return energy;
    }

    fn get_slater_orbitals (&self, basisset: &BasisSet, n: i32) -> Vec<SlaterOrbital> {

        let mut result: Vec<SlaterOrbital> = Vec::new();

        for atom in &self.atoms {

            for shell in &basisset.elements.get(&atom.atomic_number).unwrap().electron_shells {

                for orbital_coefficients in &shell.coefficients {
                    result.push(SlaterOrbital::create_slater_orbital(orbital_coefficients, n, atom.position, &shell.exponents));
                }
            }
        }

        result
    }
}

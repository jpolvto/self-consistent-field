use crate::{PI, error_function};
use crate::gaussian::Gaussian;

use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    atomic_number: i32,
    position: Vector3<f32>,
    number_of_electrons: i32,
    principal_quantum_number: i32
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Molecule {
    atoms: Vec<Atom>
}

impl Molecule {
    fn get_overlap_integral(a: &Gaussian, b: &Gaussian) -> f32 {
        let (p, difference, k, position_p) = Gaussian::gaussian_product(a, b);

        (PI / p).powf(1.5) * k
    }
    
    fn get_kinetic_integral(a: &Gaussian, b: &Gaussian) -> f32 {
        let (p, difference, k, position_p) = Gaussian::gaussian_product(a, b);
        
        (a.prefactor * b.prefactor / p) * (3.0 - 2.0 * (a.prefactor * b.prefactor / p) * difference) * (PI / p).powf(1.5) * k
    }

    fn get_potential_integral(a: &Gaussian,b: &Gaussian, atom_index: usize, molecule: Molecule) -> f32 {
        let (p, difference, k, position) = Gaussian::gaussian_product(a, b);

        (-2.0 * PI * molecule.atoms[atom_index].atomic_number as f32 / p) * k * error_function((p * (position - molecule.atoms[atom_index].position).norm_squared()).into())
    }
    
    fn get_multi_electron_integral(a: &Gaussian, b: &Gaussian, gaussian_c: &Gaussian, gaussian_d: &Gaussian) -> f32 {
        let (p, difference_ab, k_ab, position_p) = Gaussian::gaussian_product(a, b);
        let (q, difference_cd, k_cd, position_q) = Gaussian::gaussian_product(gaussian_c, gaussian_d);
    
        2.0 * PI.powf(2.5) * 1.0 / (p * q * (p + q)).sqrt() * k_ab * k_cd * error_function((p * q / (p + q) * (position_p - position_q).norm_squared()).into())
    }

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

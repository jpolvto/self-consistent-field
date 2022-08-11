use crate::basisset::BasisSet;

use nalgebra::Vector3;
use serde::{Serialize, Deserialize};
use fastapprox::faster::erf;


const PI: f32 = 3.141592654;

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

#[derive(Serialize, Deserialize, Debug)]
pub struct Gaussian {
    position: Vector3<f32>,
    prefactor: f32
}

impl Molecule {
    fn get_overlap_integral(a: &Gaussian, b: &Gaussian) -> f32 {
        let (p, difference, k, position_p) = gaussian_product(a, b);

        (PI / p).powf(1.5) * k
    }
    
    fn get_kinetic_integral(a: &Gaussian, b: &Gaussian) -> f32 {
        let (p, difference, k, position_p) = gaussian_product(a, b);
        
        (a.prefactor * b.prefactor / p) * (3.0 - 2.0 * (a.prefactor * b.prefactor / p) * difference) * (PI / p).powf(1.5) * k
    }

    fn get_potential_integral(a: &Gaussian,b: &Gaussian, atom_index: usize, molecule: Molecule) -> f32 {
        let (p, difference, k, position) = gaussian_product(a, b);

        (-2.0 * PI * molecule.atoms[atom_index].atomic_number as f32 / p) * k * error_function((p * (position - molecule.atoms[atom_index].position).norm_squared()).into())
    }
    
    fn get_multi_electron_integral(a: &Gaussian, b: &Gaussian, gaussian_c: &Gaussian, gaussian_d: &Gaussian) -> f32 {
        let (p, difference_ab, k_ab, position_p) = gaussian_product(a, b);
        let (q, difference_cd, k_cd, position_q) = gaussian_product(gaussian_c, gaussian_d);
    
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

fn gaussian_product (a: &Gaussian, b: &Gaussian) -> (f32, f32, f32, Vector3<f32>) {
    
    let p = a.prefactor + b.prefactor;
    let difference = (a.position - b.position).norm_squared();
    let k = (4.0 * a.prefactor * b.prefactor / (PI.powf(2.0))).powf(0.75) * (-a.prefactor * b.prefactor / p * difference).exp();
    let position_p = (a.prefactor * a.position + b.prefactor * b.position) / p;

    (p, difference, k, position_p)
}

fn error_function(t: f32) -> f32 {
    if t == 0.0 {
         1.0
    } else {
        (0.5 * (PI / t).powf(0.5)) * erf(t.powf(0.5))
    }
}

fn wavefunction_sto_ng(n: i32, b: BasisSet, r: f32, atomic_number: i32) -> f32 {

    let mut wavefunction: f32 = Default::default();

    for i in 0..n {
        wavefunction += b.elements.get(&0).unwrap().electron_shells[0 as usize].coefficients[0 as usize][i as usize] * ((2.0 * b.elements.get(&0).unwrap().electron_shells[0 as usize].exponents[i as usize]) / PI).powf(0.75) * (-b.elements.get(&0).unwrap().electron_shells[0 as usize].exponents[i as usize] * r.powf(2.0)).exp();
    }

    return wavefunction

}

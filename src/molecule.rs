use nalgebra::Vector3;
use serde::{Serialize, Deserialize};
use statrs::function::erf::erf;

const PI: f64 = 3.141592654;

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    atomic_number: i32,
    position: Vector3<f64>,
    number_of_electrons: i32
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Molecule {
    atoms: Vec<Atom>
}

impl Molecule {
    fn get_max_angular_momentum(&self) {
        
    }

    fn get_number_of_contracted_gaussians(&self) {
        
    }

    fn get_number_of_gaussians(&self) {
        
    }

    fn get_alpha_electrons(&self) {
        
    }

    fn get_beta_electrons(&self) {
        
    }

    fn get_nuclear_repulsion_energy(&self) -> f64 {

        let mut energy: f64 = 0.0;

        for atom1 in 0..self.atoms.len() {
            for atom2 in (atom1 + 1)..self.atoms.len() {
                energy += (self.atoms[atom1].atomic_number) as f64 / (self.atoms[atom1].position - self.atoms[atom2].position).norm();
            } 
        }
        return energy;
    }

    fn get_number_of_electrons(&self) {

    }

}

fn gaussian_product (gaussian_a:(f64, Vector3<f64>), gaussian_b:(f64, Vector3<f64>,)) -> (f64, f64, f64, Vector3<f64>) {
    let (a, position_a) = gaussian_a;
    let (b, position_b) = gaussian_b;
    
    let p = a + b;
    let difference = (position_a - position_b).norm_squared();
    let k = (4.0 * a * b / (PI.powf(2.0))).powf(0.75) * (-a * b / p * difference).exp();
    let position_p = (a * position_a + b * position_b) / p;

    (p, difference, k, position_p)
}

fn overlap_integral (orbital_a:(f64, Vector3<f64>), orbital_b:(f64, Vector3<f64>)) {
    let (p, difference, k, position_p) = gaussian_product(orbital_a, orbital_b);
    (PI / position_p).pow(1.5) * k
}

fn kinetic_integral (orbital_a:(f64, Vector3<f64>), orbital_b:(f64, Vector3<f64>)) -> f64 {
    let (p, difference, k, position_p) = gaussian_product(orbital_a, orbital_b);
    let (a, position_a) = orbital_a;
    let (b, position_b) = orbital_b;

    (a * b / p) * (3.0 - 2.0 * (a * b / p) * difference) * (PI / p).powf(1.5) * k
}

fn error_function(t: f64) -> f64 {
    if t == 0.0 {
         1.0
    } else {
        (0.5 * (PI / t).powf(0.5)) * erf(t.powf(0.5))
    }
}

fn potential_integral (gaussian_a:(f64, Vector3<f64>),gaussian_b:(f64, Vector3<f64>),atom_index: usize, molecule: Molecule) {
    let (p, difference, k, position) = gaussian_product(gaussian_a, gaussian_b);

    (-2.0 * PI * molecule.atoms[atom_index].atomic_number as f64 / position) * k * error_function((p * (position - molecule.atoms[atom_index].position).norm_squared()).into())
}

fn multi_electron_integral (gaussian_a:(f64, Vector3<f64>), gaussian_b:(f64, Vector3<f64>), gaussian_c:(f64, Vector3<f64>), gaussian_d:(f64, Vector3<f64>)) -> f64 {
    let (p, difference_ab, k_ab, position_p) = gaussian_product(gaussian_a, gaussian_b);
    let (q, difference_cd, k_cd, position_q) = gaussian_product(gaussian_c, gaussian_d);

    2.0 * PI.powf(2.5) * 1.0 / (p * q * (p + q)).sqrt() * k_ab * k_cd * error_function((p * q / (p + q) * (position_p - position_q).norm_squared()).into())
}


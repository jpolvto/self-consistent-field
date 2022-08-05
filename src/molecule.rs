use nalgebra::Vector3;
use serde::{Serialize, Deserialize};
use statrs::function::erf::erf;

const PI:f64 = 3.141592654;

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    atomic_number: i32,
    position: Vector3<i32>,
    number_of_electrons: i32
}

impl Atom {
    fn get_position(&self) -> Vector3<i32> {
        return self.position
    }

    fn get_atomic_number(&self) -> i32 {
        return self.atomic_number;
    }
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

    fn get_nuclear_repulsion_energy(&self) -> i32 {

        let mut energy:i32 = 0;

        for atom1 in 0..self.atoms.len() {
            for atom2 in (atom1 + 1)..self.atoms.len() {
                energy += self.atoms[atom1].z / (self.atoms[atom1].position - self.atoms[atom2].position).len();
            } 
        }
        return energy;
    }

    fn get_number_of_electrons(&self) {

    }

}

fn gaussian_product (gaussian_a:(i32, Vector3<i32>), gaussian_b:(i32, Vector3<i32>,)) -> (i32, i32, i32, Vector3<i32>) {
    let (a, position_a) = gaussian_a;
    let (b, position_b) = gaussian_b;
    
    let p = a + b;
    let difference = (position_a - position_b).norm_squared();
    let n = (4 * a * b / (PI.powf(2))).pow(0.75);
    let k = n * (-a * b / p * difference).exp;
    let position_p = (a * position_a + b * position_b) / p;

    (p, difference, k, position_p)
}

fn overlap_integral (orbital_a:(i32, Vector3<i32>), orbital_b:(i32, Vector3<i32>)) {
    let (p, difference, k, position_p) = gaussian_product(orbital_a, orbital_b);
    (PI / position_p).pow(1.5) * k
}

fn kinetic_integral (orbital_a:(i32, Vector3<i32>), orbital_b:(i32, Vector3<i32>)) -> i32 {
    let (p, difference, k, position_p) = gaussian_product(orbital_a, orbital_b);
    let (a, position_a) = orbital_a;
    let (b, position_b) = orbital_b;

    (a * b / p) * (3 - 2 * (a * b / p) * difference) * (PI / p).pow(1.5) * k
}

fn error_function(t: f64) -> f64 {
    if t == 0.0 {
         1.0
    } else {
        (0.5 * (PI / t).powf(0.5)) * erf(t.powf(0.5))
    }
}

fn potential_integral (gaussian_a:(i32, Vector3<i32>),gaussian_b:(i32, Vector3<i32>),atom_index: usize, molecule: Molecule) {
    let (p, difference, k, position) = gaussian_product(gaussian_a, gaussian_b);
    let atom_position = molecule.atoms[atom_index].position;
    let atomic_number = molecule.atoms[atom_index].atomic_number;

    (-2 * PI * atomic_number / position) * k * error_function(p * (position - atom_position).norm_squared())
}

fn multi_electron_integral (gaussian_a:(i32, Vector3<i32>), gaussian_b:(i32, Vector3<i32>), gaussian_c:(i32, Vector3<i32>), gaussian_d:(i32, Vector3<i32>)) -> Vector3<i32> {
    let (p, difference_ab, k_ab, position_p) = gaussian_product(gaussian_a, gaussian_b);
    let (q, difference_cd, k_cd, position_q) = gaussian_product(gaussian_c, gaussian_d);

    2 * PI.powf(2.5) * (p * q * (p + q).pow(0.5).pow(-1)) * k_ab * k_cd * error_function(p * q / (p + q) * (position_p - position_q).norm_squared())
}


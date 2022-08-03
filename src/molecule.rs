use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    z: i32,
    position: Vector3<i32>,
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

fn gaussian_product (gaussian_a:(i32, Vector3<i32>), gaussian_b:(i32, Vector3<i32>,)) {
    let a = gaussian_a.0;
    let Ra = gaussian_a.1;
    let b = gaussian_b.0;
    let Rb = gaussian_b.1;

    let pi = 3.141592653;

    let p = a + b;

    let difference = (Ra - Rb).norm_squared();

    let N = (4*a*b/(pi.pow(2))).pow(0.75);
    let K = N*exp(-a*b/p*difference)
    let Rp = (a*Ra + b*Rb)/p

    return p, diff, K, Rp;
}


fn overlap_integral (orbital_a, orbital b) {
    p, difference, K, Rp = gaussian_product(orbital_a, orbital_b)
   let prefactor = (pi/p).pow(1.5);
   prefactor * K

}   

fn kinetic_integral (orbital_a, orbital b) {
    p, difference, K, Rp = gaussian_product(orbital_a, orbital_b)
    let prefactor = (pi/p).pow(1.5);

    let a = gaussian_a.0;
    let Ra = gaussian_a.1;
    let b = gaussian_b.0;
    let Rb = gaussian_b.1;

    reduced_exponent = a*b/p
    reduced_exponent*(3-2*reduced_exponent*difference)*prefactor*K
}

fn error_function(t) {
    if t == 0 {
         1
    } else {
        (0.5*(pi/t).pow(0.5))*erf(t.pow(0.5))
    }
}

fn potential_integral
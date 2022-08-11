use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::basisset::BasisSet;
use crate::PI;

#[derive(Serialize, Deserialize, Debug)]
pub struct Gaussian {
    pub position: Vector3<f32>,
    pub prefactor: f32
}

impl Gaussian {
    fn wavefunction_sto_ng(n: i32, b: BasisSet, r: f32, atomic_number: i32) -> f32 {

        let mut wavefunction: f32 = Default::default();
    
        for i in 0..n {
            wavefunction += b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].coefficients[usize::MIN][i as usize] * ((2.0 * b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].exponents[i as usize]) / PI).powf(0.75) * (-b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].exponents[i as usize] * r.powf(2.0)).exp();
        }
    
        return wavefunction
    
    }

    pub fn gaussian_product (a: &Gaussian, b: &Gaussian) -> (f32, f32, f32, Vector3<f32>) {
    
        let p = a.prefactor + b.prefactor;
        let difference = (a.position - b.position).norm_squared();
        let k = (4.0 * a.prefactor * b.prefactor / (PI.powf(2.0))).powf(0.75) * (-a.prefactor * b.prefactor / p * difference).exp();
        let position_p = (a.prefactor * a.position + b.prefactor * b.position) / p;
    
        (p, difference, k, position_p)
    }
}


use std::f32::consts::PI;
use fastapprox::faster::erf;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::{basisset::BasisSet, molecule::Atom};

#[derive(Serialize, Deserialize, Debug)]
pub struct Gaussian {
    pub center: Vector3<f32>,
    pub exponent: f32
}

impl Gaussian {
    fn wavefunction_sto_ng(n: i32, b: BasisSet, r: f32, atomic_number: i32) -> f32 {

        let mut wavefunction: f32 = Default::default();
    
        for i in 0..n {
            wavefunction += b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].coefficients[usize::MIN][i as usize] * ((2.0 * b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].exponents[i as usize]) / PI).powf(0.75) * (-b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].exponents[i as usize] * r.powf(2.0)).exp();
        }
        return wavefunction
    }

    pub fn gaussian_product (a: &Gaussian, b: &Gaussian) -> (Gaussian, f32, f32) {

        let difference = (a.center - b.center).norm_squared();
        let exponent = a.exponent + b.exponent;
        let k = (4.0 * a.exponent * b.exponent / (PI.powf(2.0))).powf(0.75) * (-a.exponent * b.exponent / exponent * difference).exp();
        let center = (a.exponent * a.center + b.exponent * b.center) / exponent;

        (Gaussian {center, exponent}, difference, k)
    }

    fn get_overlap_integral(a: &Gaussian, b: &Gaussian) -> f32 {
        let (c, difference, k) = Gaussian::gaussian_product(a, b);

        (PI / c.exponent).powf(1.5) * k
    }
    
    fn get_kinetic_integral(a: &Gaussian, b: &Gaussian) -> f32 {
        let (c, difference, k) = Gaussian::gaussian_product(a, b);
        
        (a.exponent * b.exponent / c.exponent) * (3.0 - 2.0 * (a.exponent * b.exponent / c.exponent) * difference) * (PI / c.exponent).powf(1.5) * k
    }

    fn get_potential_integral(a: &Gaussian,b: &Gaussian, atom: Atom) -> f32 {
        let (c, difference, k) = Gaussian::gaussian_product(a, b);

        (-2.0 * PI * atom.atomic_number as f32 / c.exponent) * k * Gaussian::error_function((c.exponent * (c.center - atom.position).norm_squared()).into())
    }
    
    fn get_multi_electron_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, d: &Gaussian) -> f32 {
        let (p, difference_ab, k_ab) = Gaussian::gaussian_product(a, b);
        let (q, difference_cd, k_cd) = Gaussian::gaussian_product(c, d);
    
        2.0 * PI.powf(2.5) * 1.0 / (p.exponent * q.exponent * (p.exponent + q.exponent)).sqrt() * k_ab * k_cd * Gaussian::error_function((p.exponent * q.exponent / (p.exponent + q.exponent) * (p.center - q.center).norm_squared()).into())
    }

    fn error_function(t: f32) -> f32 {
        if t == 0.0 {
             1.0
        } else {
            (0.5 * (PI / t).powf(0.5)) * erf(t.powf(0.5))
        }
    }
    
}

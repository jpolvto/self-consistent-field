use std::f32::consts::PI;

use fastapprox::fast::erf;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::molecule::Atom;

#[derive(Debug, Serialize, Deserialize)]
pub struct Gaussian {
    //nuclear coordinates
    pub center: Vector3<f32>,

    //Gaussian orbital coefficient
    pub coefficient: f32,
}

impl Gaussian {

    pub fn gaussian_product(&self, other: &Gaussian) -> (Gaussian, f32, f32) {

        // function to calculate a gaussian product, and to have some convenient variables

        let diff = (self.center-other.center).norm_squared();
        let coefficient = self.coefficient + other.coefficient;
        let n = (4.0*self.coefficient*other.coefficient/(PI.powf(2.0))).powf(0.75);
        let k = n*(-self.coefficient*other.coefficient/coefficient*diff).exp();

        let center = ((self.coefficient*self.center + other.coefficient*other.center))/coefficient;

        (Gaussian {center, coefficient}, diff, k)
    }

    pub fn overlap_integral(&self, other: &Gaussian) -> f32 {

        //Calculates the overlap integral between two gaussian functions

        let (p, _diff_ab, k) = Gaussian::gaussian_product(self, other);

        (PI/p.coefficient).powf(1.5)*k
    
        //normalization constant

    }
    
    pub fn kinetic_energy_integral(&self, other: &Gaussian) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        let (p, diff_ab, k) = Gaussian::gaussian_product(self, other);
        let reduced_exponent = self.coefficient*other.coefficient/p.coefficient;

        reduced_exponent*(3.0-2.0*reduced_exponent*diff_ab)*(PI/p.coefficient).powf(1.5)*k

    }

    pub fn nuclear_attraction_integral(&self, other: &Gaussian, atom: &Atom) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        let (p, _diff_ab, k) = Gaussian::gaussian_product(self, other);
        let mut matrix_element = (-2.0*PI*atom.atomic_number as f32/p.coefficient)*k;

        if let Some(i) = Gaussian::f0(p.coefficient*(p.center-atom.position).norm_squared()) {
            matrix_element *= i;
        }
        matrix_element
    }
    
    pub fn two_electron_integral(&self, a: &Gaussian, b: &Gaussian, c: &Gaussian) -> f32 {

        /*
        Calculate two electron integrals
        alpha, beta, delta, gamma are the coefficients of the gaussian orbitals
        r_p is the distance between a and b, r_q is the distance between c and d
        */

        let (p, _diff_ab, k_ab) = self.gaussian_product(a);
        let (q, _diff_cd, k_cd) = b.gaussian_product(c);

        let mut matrix_element = k_ab*k_cd*2.0*PI.powf(2.5)/(p.coefficient*q.coefficient*(p.coefficient+q.coefficient).sqrt());

        if let Some(i) = Gaussian::f0(p.coefficient*q.coefficient/(p.coefficient+q.coefficient)*(p.center-q.center).norm_squared()) {
            matrix_element *= i;
        }
        matrix_element
    }

    fn f0(t: f32) -> Option<f32> {

        // single electron error function
        
        if t.abs() < 1e-8 {
             return None
        } else {
            return Some(0.5*(PI/t).sqrt()*erf(t.sqrt()))
        }
    }    
}
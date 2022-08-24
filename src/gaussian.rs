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

    pub fn gaussian_product(a: &Gaussian, b: &Gaussian) -> (Gaussian, f32, f32) {

        // function to calculate a gaussian product, and to have some convenient variables

        let r_ab = (a.center-b.center).norm_squared();
        let coefficient = a.coefficient+b.coefficient;
        let k = (4.0*a.coefficient*b.coefficient/(PI.powf(2.0))).powf(0.75)*(-a.coefficient*b.coefficient/coefficient*r_ab).exp();
        let center = (a.coefficient*a.center+b.coefficient*b.center)/coefficient;

        (Gaussian {center, coefficient}, r_ab, k)
    }

    pub fn overlap_integral(a: &Gaussian, b: &Gaussian) -> f32 {

        //Calculates the overlap integral between two gaussian functions

        let (p, _r_ab, k) = Gaussian::gaussian_product(a, b);

        (PI/p.coefficient).powf(1.5)*k
    
        //normalization constant

    }
    
    pub fn kinetic_energy_integral(a: &Gaussian, b: &Gaussian) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        let (p, r_ab, k) = Gaussian::gaussian_product(a, b);

        (a.coefficient*b.coefficient/p.coefficient)*(3.0-2.0*
        (a.coefficient*b.coefficient/p.coefficient)*r_ab)*(PI/p.coefficient).powf(1.5)*k

    }

    pub fn nuclear_attraction_integral(a: &Gaussian, b: &Gaussian, atom: &Atom) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        let (p, _r_ab, k) = Gaussian::gaussian_product(a, b);

        let mut matrix_element = (-2.0*PI*atom.atomic_number as f32/p.coefficient)*k;

        if let Some(i) = Gaussian::f0(p.coefficient*(p.center-atom.position).norm_squared()) {
            matrix_element *= i;
        }
        matrix_element
    }
    
    pub fn two_electron_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, d: &Gaussian) -> f32 {

        /*
        Calculate two electron integrals
        alpha, beta, delta, gamma are the coefficients of the gaussian orbitals
        r_p is the distance between a and b, r_q is the distance between c and d
        */

        let (p, _r_ab, k_ab) = Gaussian::gaussian_product(a, b);
        let (q, _r_cd, k_cd) = Gaussian::gaussian_product(c, d);

        let mut matrix_element = 2.0*PI.powf(2.5)/(p.coefficient*q.coefficient*(p.coefficient+q.coefficient)).sqrt()*k_ab*k_cd;

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
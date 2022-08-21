use std::f32::consts::PI;
use fastapprox::faster::erf;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::molecule::Atom;

#[derive(Debug, Serialize, Deserialize)]
pub struct Gaussian {
    //nuclear coordinates
    pub center: Vector3<f32>,

    //Gaussian orbital coefficient
    pub coefficient: f32,

    //Slater orbital exponent
    pub exponent: f32
}

impl Gaussian {

    pub fn gaussian_product(a: &Gaussian, b: &Gaussian) -> (Gaussian, f32, f32) {

        // function to calculate a gaussian product, and to have some convenient variables

        let r_ab = (a.center-b.center).norm_squared();
        let coefficient = a.coefficient+b.coefficient;
        let k = (4.0*a.coefficient*b.coefficient/(PI.powf(2.0))).powf(0.75)*(-a.coefficient*b.coefficient/coefficient*r_ab).exp();
        let center = (a.coefficient*a.center+b.coefficient*b.center)/coefficient;
        let exponent = a.exponent*b.exponent;

        (Gaussian {center, coefficient, exponent}, r_ab, k)
    }

    pub fn overlap_integral(p: &Gaussian, k: f32) -> f32 {

        //Calculates the overlap between two gaussian functions

        (PI/p.coefficient).powf(1.5)*k*p.exponent

    }
    
    pub fn kinetic_energy_integral(a: &Gaussian, b: &Gaussian, p: &Gaussian, r_ab: f32, k: f32) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        (a.coefficient*b.coefficient/p.coefficient)*(3.0-2.0*
        (a.coefficient*b.coefficient/p.coefficient)*r_ab)*(PI/p.coefficient).powf(1.5)*k*p.exponent

    }

    pub fn nuclear_attraction_integral(atom: &Atom, p: &Gaussian, k: f32) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        (-2.0*PI*atom.atomic_number as f32/p.coefficient)*k*
        Gaussian::f_zero(p.coefficient*(p.center-atom.position).norm_squared())*p.exponent

    }
    
    pub fn two_electron_integral(p: &Gaussian, q: &Gaussian, k_ab: f32, k_cd: f32) -> f32 {

        /*
        Calculate two electron integrals
        a.coefficient,b.coefficient,c.coefficient,d.coefficient are the coefficients alpha, beta, etc.
        rab2 equals squared distance between centre a.coefficient and centre b.coefficient
        */

        2.0*PI.powf(2.5)/(p.coefficient*q.coefficient*(p.coefficient+q.coefficient)).sqrt()*k_ab*k_cd*
        Gaussian::f_zero(p.coefficient*q.coefficient/(p.coefficient+q.coefficient)*
        (p.center-q.center).norm_squared())*p.exponent*q.exponent

    }

    fn f_zero(t: f32) -> f32 {

        // single elctron error function
        
        if t == 0.0 {
             1.0
        } else {
            (0.5*(PI/t).powf(0.5))*erf(t.powf(0.5))
        }
    }
    
}


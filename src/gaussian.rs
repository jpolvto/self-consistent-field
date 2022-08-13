use std::f32::consts::PI;
use fastapprox::faster::erf;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::{basisset::BasisSet, molecule::{Atom, Molecule}};

#[derive(Debug)]
pub struct Gaussian {
    //nuclear coordinates
    pub center: Vector3<f32>,

    //Gaussian orbital exponent
    pub exponent: f32
}

//Slater Type Orbital fit with N primative gausians (STO-NG) type basis
#[derive(Debug)]
pub struct SlaterOrbital {
    //linear combination of n Gaussian orbitals
    n: i32,

    //contraction coefficients
    zetas: Vec<f32>,

    //primitive Gaussians
    gaussians: Vec<Gaussian>
}

impl Gaussian {

    fn gaussian_product(a: &Gaussian, b: &Gaussian) -> (Gaussian, f32, f32) {

        // function to calculate a gaussian product, and to have some convenient variables

        let r_ab = (a.center-b.center).norm_squared();
        let exponent = a.exponent+b.exponent;
        let k = (4.0*a.exponent*b.exponent/(PI.powf(2.0))).powf(0.75)*(-a.exponent*b.exponent/exponent*r_ab).exp();
        let center = (a.exponent*a.center+b.exponent*b.center)/exponent;

        (Gaussian {center, exponent}, r_ab, k)
    }

    fn get_overlap_integral(c: &Gaussian, k: f32) -> f32 {

        //Calculates the overlap between two gaussian functions

        (PI/c.exponent).powf(1.5)*k

    }
    
    fn get_kinetic_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, r_ab: f32, k: f32) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        (a.exponent*b.exponent/c.exponent)*(3.0-2.0*
        (a.exponent*b.exponent/c.exponent)*r_ab)*(PI/c.exponent).powf(1.5)*k

    }

    fn get_potential_integral(atom: Atom, c: &Gaussian, k: f32) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        (-2.0*PI*atom.atomic_number as f32/c.exponent)*k*
        Gaussian::f_zero(c.exponent*(c.center-atom.position).norm_squared())

    }
    
    fn get_multi_electron_integral(p: &Gaussian, q: &Gaussian, k_ab: f32, k_cd: f32) -> f32 {

        /*
        Calculate two electron integrals
        a.exponent,b.exponent,c.exponent,d.exponent are the exponents alpha, beta, etc.
        rab2 equals squared distance between centre a.exponent and centre b.exponent
        */

        2.0*PI.powf(2.5)/(p.exponent*q.exponent*(p.exponent+q.exponent)).sqrt()*k_ab*k_cd*
        Gaussian::f_zero(p.exponent*q.exponent/(p.exponent+q.exponent)*
        (p.center-q.center).norm_squared())

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


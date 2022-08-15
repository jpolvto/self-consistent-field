use std::f32::consts::PI;
use fastapprox::faster::erf;
use nalgebra::Vector3;

use crate::molecule::Atom;

//Slater Type Orbital fit with N primative gausians (STO-NG) type basis
#[derive(Debug)]
pub struct SlaterOrbital {
    //linear combination of n Gaussian orbitals
    n: i32,

    //contraction coefficients
    exponents: Vec<f32>,

    //primitive Gaussians
    gaussians: Vec<Gaussian>
}

impl SlaterOrbital {

    pub fn create_slater_orbital(orbital_coefficients: &Vec<f32>, n: i32, center: Vector3<f32>, exponents: &Vec<f32>) -> SlaterOrbital {

        let mut gaussians: Vec<Gaussian> = Vec::new();
            
        for coefficient in orbital_coefficients {
            gaussians.push( Gaussian { center, coefficient: coefficient.clone() })
        }

        SlaterOrbital { n, exponents: exponents.clone(), gaussians }
    }
}

#[derive(Debug)]
pub struct Gaussian {
    //nuclear coordinates
    pub center: Vector3<f32>,

    //Gaussian orbital coefficient
    pub coefficient: f32
}

impl Gaussian {

    fn gaussian_product(a: &Gaussian, b: &Gaussian) -> (Gaussian, f32, f32) {

        // function to calculate a gaussian product, and to have some convenient variables

        let r_ab = (a.center-b.center).norm_squared();
        let coefficient = a.coefficient+b.coefficient;
        let k = (4.0*a.coefficient*b.coefficient/(PI.powf(2.0))).powf(0.75)*(-a.coefficient*b.coefficient/coefficient*r_ab).exp();
        let center = (a.coefficient*a.center+b.coefficient*b.center)/coefficient;

        (Gaussian {center, coefficient}, r_ab, k)
    }

    fn get_overlap_integral(c: &Gaussian, k: f32) -> f32 {

        //Calculates the overlap between two gaussian functions

        (PI/c.coefficient).powf(1.5)*k

    }
    
    fn get_kinetic_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, r_ab: f32, k: f32) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        (a.coefficient*b.coefficient/c.coefficient)*(3.0-2.0*
        (a.coefficient*b.coefficient/c.coefficient)*r_ab)*(PI/c.coefficient).powf(1.5)*k

    }

    fn get_potential_integral(atom: &Atom, c: &Gaussian, k: f32) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        (-2.0*PI*atom.atomic_number as f32/c.coefficient)*k*
        Gaussian::f_zero(c.coefficient*(c.center-atom.position).norm_squared())

    }
    
    fn get_multi_electron_integral(p: &Gaussian, q: &Gaussian, k_ab: f32, k_cd: f32) -> f32 {

        /*
        Calculate two electron integrals
        a.coefficient,b.coefficient,c.coefficient,d.coefficient are the coefficients alpha, beta, etc.
        rab2 equals squared distance between centre a.coefficient and centre b.coefficient
        */

        2.0*PI.powf(2.5)/(p.coefficient*q.coefficient*(p.coefficient+q.coefficient)).sqrt()*k_ab*k_cd*
        Gaussian::f_zero(p.coefficient*q.coefficient/(p.coefficient+q.coefficient)*
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


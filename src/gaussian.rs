use std::f32::consts::PI;
use fastapprox::faster::erf;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::{basisset::BasisSet, molecule::{Atom, Molecule}};

#[derive(Serialize, Deserialize, Debug)]
pub struct Gaussian {
    pub center: Vector3<f32>,
    pub exponent: f32
}

impl Gaussian {
    fn wavefunction_sto_ng(n: i32, b: BasisSet, r: f32, atomic_number: i32) -> f32 {

        let mut wavefunction: f32 = Default::default();
    
        for i in 0..n {
            wavefunction += b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].coefficients[usize::MIN][i as usize]*
            ((2.0*b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].exponents[i as usize])/PI).powf(0.75)*
            (-b.elements.get(&atomic_number).unwrap().electron_shells[usize::MIN].exponents[i as usize]*r.powf(2.0)).exp();
        }
        return wavefunction
    }

    fn gaussian_product(a: &Gaussian, b: &Gaussian) -> (Gaussian, f32, f32) {

        let difference = (a.center-b.center).norm_squared();
        let exponent = a.exponent+b.exponent;
        let k = (4.0*a.exponent*b.exponent/(PI.powf(2.0))).powf(0.75)*(-a.exponent*b.exponent/exponent*difference).exp();
        let center = (a.exponent*a.center+b.exponent*b.center)/exponent;

        (Gaussian {center, exponent}, difference, k)
    }

    fn get_overlap_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, k: f32) -> f32 {

        //Calculates the overlap between two gaussian functions

        (PI/c.exponent).powf(1.5)*k

    }
    
    fn get_kinetic_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, difference: f32, k: f32) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        (a.exponent*b.exponent/c.exponent)*(3.0-2.0*
        (a.exponent*b.exponent/c.exponent)*difference)*(PI/c.exponent).powf(1.5)*k

    }

    fn get_potential_integral(a: &Gaussian,b: &Gaussian, atom: Atom, c: &Gaussian, k: f32) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        (-2.0*PI*atom.atomic_number as f32/c.exponent)*k*
        Gaussian::error_function(c.exponent*(c.center-atom.position).norm_squared())

    }
    
    fn get_multi_electron_integral(p: &Gaussian, q: &Gaussian, kab: f32, kcd: f32) -> f32 {

        /*
        Calculate two electron integrals
        a.exponent,b.exponent,c.exponent,d.exponent are the exponents alpha, beta, etc.
        rab2 equals squared distance between centre a.exponent and centre b.exponent
        */

        2.0*PI.powf(2.5)/(p.exponent*q.exponent*(p.exponent+q.exponent)).sqrt()*kab*kcd*
        Gaussian::error_function(p.exponent*q.exponent/(p.exponent+q.exponent)*
        (p.center-q.center).norm_squared())

    }

    fn error_function(t: f32) -> f32 {
        if t == 0.0 {
             1.0
        } else {
            (0.5*(PI/t).powf(0.5))*erf(t.powf(0.5))
        }
    }
    
}


fn overlap_integral(a: Gaussian, b: Gaussian, rab2: f32) -> f32 {

    //Calculates the overlap between two gaussian functions 
    
    (PI/(a.exponent+b.exponent)).powf(1.5)*(-a.exponent*b.exponent*rab2/(a.exponent+b.exponent)).exp()
}


fn kinetic_integral(a: Gaussian, b: Gaussian, rab2: f32) -> f32 {

    //Calculates the kinetic energy integrals for un-normalised primitives
    
    a.exponent*b.exponent/(a.exponent+b.exponent)*(3.0-2.0*a.exponent*b.exponent*rab2/(a.exponent+b.exponent))*
    (PI/(a.exponent+b.exponent)).powf(1.5)*(-a.exponent*b.exponent*rab2/(a.exponent+b.exponent)).exp()
}

 fn nuclear_attraction_integral(a: Gaussian, b: Gaussian, rab2: f32, rcp2: f32, atomic_number: i32) -> f32 {

    //Calculates the un-normalised nuclear attraction integrals

    -2.0*PI/(a.exponent+b.exponent)*Gaussian::error_function((a.exponent+b.exponent)*rcp2)*
    (-a.exponent*b.exponent*rab2/(a.exponent+b.exponent)).exp()*atomic_number as f32
 }

fn two_electron_integral(a: Gaussian, b: Gaussian, c: Gaussian, d: Gaussian, rab2: f32, rcd2: f32, rpq2: f32) -> f32 {

    /*
    Calculate two electron integrals
    a.exponent,b.exponent,c.exponent,d.exponent are the exponents alpha, beta, etc.
    rab2 equals squared distance between centre a.exponent and centre b.exponent
    */

    2.0*(PI.powf(2.5)/((a.exponent+b.exponent)*(c.exponent+d.exponent)*(a.exponent+b.exponent+c.exponent+d.exponent).sqrt())*
    Gaussian::error_function((a.exponent+b.exponent)*(c.exponent+d.exponent)*rpq2/(a.exponent+b.exponent+c.exponent+d.exponent))*
    (-a.exponent*b.exponent*rab2/(a.exponent+b.exponent)).exp()-c.exponent*d.exponent*rcd2/(c.exponent+d.exponent))
}


fn create_gaussian_orbitals(molecule: Molecule, basisset: BasisSet) {
    for atom in molecule.atoms {

    }
}

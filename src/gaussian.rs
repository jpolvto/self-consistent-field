use std::f32::consts::PI;
use fastapprox::faster::erf;
use nalgebra::Vector3;
use serde::{Serialize, Deserialize};

use crate::molecule::Atom;

//Slater Type Orbital fit with N primative gausians (STO-NG) type basis
#[derive(Debug, Serialize, Deserialize)]
pub struct Orbital {
    //contraction coefficients
    exponents: Vec<f32>,

    //primitive Gaussians
    gaussians: Vec<Gaussian>
}

impl Orbital {

    pub fn create_orbital(orbital_coefficients: &Vec<f32>, center: Vector3<f32>, exponents: &Vec<f32>) -> Orbital {

        let mut gaussians: Vec<Gaussian> = Vec::new();
        for coefficient in orbital_coefficients {
            gaussians.push( Gaussian { center, coefficient: coefficient.clone() })
        }
        Orbital { exponents: exponents.clone(), gaussians }
    }

    pub fn kinetic_energy(a: &Orbital, b: &Orbital) -> f32{

        let mut kinetic_energy: f32 = Default::default();

        for i in 0..2 {
            for j in 0..2 {
                let gaussian_a = a.gaussians.get(i).unwrap();
                let gaussian_b = b.gaussians.get(j).unwrap();

                let exp_a = a.exponents.get(i).unwrap();
                let exp_b = b.exponents.get(j).unwrap();
                let (p, r_ab, k_ab) = Gaussian::gaussian_product(gaussian_a, gaussian_b);
                kinetic_energy += Gaussian::kinetic_energy_integral(&gaussian_a, &gaussian_b, &p, r_ab, k_ab)*exp_a*exp_b;
            }
        }
        kinetic_energy
    }

    pub fn overlap_energy(a: &Orbital, b: &Orbital) -> f32{
        let mut overlap_energy: f32 = Default::default();

        for i in 0..2 {
            for j in 0..2 {
                let gaussian_a = a.gaussians.get(i).unwrap();
                let gaussian_b = b.gaussians.get(j).unwrap();

                let exp_a = a.exponents.get(i).unwrap();
                let exp_b = b.exponents.get(j).unwrap();

                let (p, _r_ab, k_ab) = Gaussian::gaussian_product(&gaussian_a, &gaussian_b);
                overlap_energy += Gaussian::overlap_integral(&p, k_ab)*exp_a*exp_b;
            }
        }
        overlap_energy
    }

    pub fn nuclear_repulsion_energy(a: &Orbital, b: &Orbital, atom: &Atom) -> f32 {
        let mut nuclear_repulsion_energy: f32 = Default::default();

        for i in 0..2 {
            for j in 0..2 {
                let gaussian_a = a.gaussians.get(i).unwrap();
                let gaussian_b = b.gaussians.get(j).unwrap();

                let exp_a = a.exponents.get(i).unwrap();
                let exp_b = b.exponents.get(j).unwrap();

                let (p, _r_ab, k_ab) = Gaussian::gaussian_product(&gaussian_a, &gaussian_b);
                nuclear_repulsion_energy += Gaussian::nuclear_attraction_integral(&atom, &p, k_ab)*exp_a*exp_b;
            }
        }
        nuclear_repulsion_energy
    }

    pub fn two_electron_energy(a: &Orbital, b: &Orbital, c: &Orbital, d: &Orbital) -> f32 {
        let mut two_electron_energy: f32 = Default::default();

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        let gaussian_a = a.gaussians.get(i).unwrap();
                        let gaussian_b = b.gaussians.get(j).unwrap();
                        let gaussian_c = c.gaussians.get(k).unwrap();
                        let gaussian_d = d.gaussians.get(l).unwrap();

                        let exp_a = a.exponents.get(i).unwrap();
                        let exp_b = b.exponents.get(j).unwrap();
                        let exp_c = c.exponents.get(k).unwrap();
                        let exp_d = d.exponents.get(l).unwrap();

                        let (p, _r_ab, k_ab) = Gaussian::gaussian_product(&gaussian_a, &gaussian_b);
                        let (q, _r_cd, k_cd) = Gaussian::gaussian_product(&gaussian_c, &gaussian_d);    
                        two_electron_energy += Gaussian::two_electron_integral(&p, &q, k_ab, k_cd)*exp_a*exp_b*exp_c*exp_d;
                    }
                }
            }
        }
        two_electron_energy
    }

}

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

    pub fn overlap_integral(p: &Gaussian, k: f32) -> f32 {

        //Calculates the overlap between two gaussian functions

        (PI/p.coefficient).powf(1.5)*k

    }
    
    pub fn kinetic_energy_integral(a: &Gaussian, b: &Gaussian, p: &Gaussian, r_ab: f32, k: f32) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        (a.coefficient*b.coefficient/p.coefficient)*(3.0-2.0*
        (a.coefficient*b.coefficient/p.coefficient)*r_ab)*(PI/p.coefficient).powf(1.5)*k

    }

    pub fn nuclear_attraction_integral(atom: &Atom, p: &Gaussian, k: f32) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        (-2.0*PI*atom.atomic_number as f32/p.coefficient)*k*
        Gaussian::f_zero(p.coefficient*(p.center-atom.position).norm_squared())

    }
    
    pub fn two_electron_integral(p: &Gaussian, q: &Gaussian, k_ab: f32, k_cd: f32) -> f32 {

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

        // single electron error function
        
        if t == 0.0 {
             1.0
        } else {
            (0.5*(PI/t).powf(0.5))*erf(t.powf(0.5))
        }
    }
    
}


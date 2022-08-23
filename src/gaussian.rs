use std::{f32::consts::PI};
use fastapprox::faster::erf;
use nalgebra::{Vector3};
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

    pub fn two_center_contraction<F>(a: &Orbital, b: &Orbital, atom: &Atom, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Atom) -> f32 {
        let mut total: f32 = Default::default();

        for i in 0..2 {
            for j in 0..2 {
                let gaussian_a = a.gaussians.get(i).unwrap();
                let gaussian_b = b.gaussians.get(j).unwrap();

                let exp_a = a.exponents.get(i).unwrap();
                let exp_b = b.exponents.get(j).unwrap();

                total += exp_a*exp_b*integral(gaussian_a, gaussian_b, atom)
            }
        }
        total
    }

    pub fn four_center_contraction<F>(a: &Orbital, b: &Orbital, c: &Orbital, d: &Orbital, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Gaussian, &Gaussian) -> f32 {
        let mut total: f32 = Default::default();
        
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

                        total += exp_a*exp_b*exp_c*exp_d*integral(gaussian_a, gaussian_b, gaussian_c, gaussian_d)
                    }
                }
            }
        }
        total
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

    pub fn overlap_integral(a: &Gaussian, b: &Gaussian, _atom: &Atom) -> f32 {

        //Calculates the overlap integral between two gaussian functions

        let alpha = a.coefficient;
        let beta = b.coefficient;
        let r_a = a.center;
        let r_b = b.center;
    
        //normalization constant
        let n = (2.0*alpha/PI).powf(0.75)*(2.0*beta/PI).powf(0.75);
    
        let mut matrix_element= n*(PI/(alpha+beta)).powf(1.5);
        matrix_element *= -alpha*beta/(alpha+beta).exp() * (r_a-r_b).norm_squared();
        matrix_element

    }
    
    pub fn kinetic_energy_integral(a: &Gaussian, b: &Gaussian, _atom: &Atom) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        let alpha = a.coefficient;
        let beta = b.coefficient;
        let r_a = a.center;
        let r_b = b.center;
    
        //normalization constant
        let n = (2.0*alpha/PI).powf(0.75)*(2.0*beta/PI).powf(0.75);
    
        let mut matrix_element = n*alpha*beta/(alpha+beta);
        matrix_element *= 3.0-2.0*alpha*beta/((alpha+beta)/(r_a-r_b).norm_squared());
        matrix_element *= (PI/(alpha+beta)).powf(1.5);
        matrix_element *= (-alpha*beta/(alpha+beta)*(r_a-r_b).norm_squared()).exp();

        matrix_element

    }

    pub fn nuclear_attraction_integral(a: &Gaussian, b: &Gaussian, atom: &Atom) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals

        let alpha = a.coefficient;
        let beta = b.coefficient;
        let r_a = a.center;
        let r_b = b.center;
        let r_p = (alpha*r_a + beta*r_b)/(alpha + beta);
    
        //normalization constant
        let n = (2.0*alpha/PI).powf(0.75)*(2.0*beta/PI).powf(0.75);

        let mut matrix_element = n*-2.0*PI/(alpha+beta)*atom.atomic_number as f32;
        matrix_element *= (-alpha*beta/(alpha+beta)*(r_a-r_b).norm_squared()).exp();
    
        let t = (alpha+beta)*(r_p-atom.position).norm_squared();

        match Self::f_zero(t) {
            Some(i) => {
                return matrix_element*i
            }
            None => {
                return matrix_element
            }
        }

    }
    
    pub fn two_electron_integral(a: &Gaussian, b: &Gaussian, c: &Gaussian, d: &Gaussian) -> f32 {

        /*
        Calculate two electron integrals
        alpha, beta, delta, gamma are the coefficients of the gaussian orbitals
        r_p is the distance between a and b, r_q is the distance between c and d
        */

        let alpha = a.coefficient;
        let beta = b.coefficient;
        let gamma = c.coefficient;
        let delta = d.coefficient;

        let r_a = a.center;
        let r_b = b.center;
        let r_c = c.center;
        let r_d = d.center;
        let r_p = (alpha*r_a + beta*r_b)/(alpha + beta);
        let r_q = (gamma*r_c + delta*r_d)/(gamma + delta);

        //normalization constant

        let mut n = (2.0*alpha/PI).powf(0.75)*(2.0*beta/PI).powf(0.75);
        n *= (2.0*gamma/PI).powf(0.75)*(2.0*delta/PI).powf(0.75);
    
        let mut matrix_element = n*2.0*PI.powf(2.6);
        matrix_element /= (alpha+beta)*(gamma+delta)*(alpha+beta+gamma+delta).sqrt();
        matrix_element *= (-alpha*beta/(alpha+beta)*(r_a-r_b).norm_squared()-gamma*delta/(gamma+delta)*(r_c-r_d).norm_squared()).exp();

        let t = (alpha+beta)*(gamma+delta)/(alpha+beta+gamma+delta)*(r_p-r_q).norm_squared();

        match Self::f_zero(t) {
            Some(i) => {
                return matrix_element*i
            }
            None => {
                return matrix_element
            }
        }

    }

    fn f_zero(t: f32) -> Option<f32> {

        // single electron error function
        
        if t.abs() < 1e-8 {
             return None
        } else {
            return Some(0.5*(PI/t).sqrt()*erf(t.sqrt()))
        }
    }    
}
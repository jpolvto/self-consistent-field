use std::{f32::consts::PI};
use fastapprox::fast::erf;
use serde::{Serialize, Deserialize};

use crate::molecule::Atom;

#[derive(Debug, Serialize, Deserialize)]
pub struct Gaussian {
    //nuclear coordinates
    pub center: f32,

    //Gaussian orbital coefficient
    pub coefficient: f32,
}

impl Gaussian {

    pub fn overlap_integral(&self, other: &Gaussian) -> f32 {

        //Calculates the overlap integral between two gaussian functions

        // normalization constant
        let n = (2.0 *self.coefficient/PI).powf(0.75) * (2.0 *other.coefficient/PI).powf(0.75);

        let mut matrix_element = n * (PI/(self.coefficient+other.coefficient)).powf(0.75); 
        matrix_element *= -self.coefficient*other.coefficient/(self.coefficient+other.coefficient).exp() * (self.center-other.center).abs().powf(2.0);
        matrix_element
    }
    
    pub fn kinetic_energy_integral(&self, other: &Gaussian) -> f32 {

        //Calculates the kinetic energy integrals for un-normalised primitives

        // normalization constant
        let n = (2.0 *self.coefficient/PI).powf(0.75) * (2.0 *other.coefficient/PI).powf(0.75);

        let mut matrix_element = n * self.coefficient*other.coefficient/(self.coefficient+other.coefficient);
        matrix_element *= 3.0-2.0*self.coefficient*other.coefficient/((self.coefficient+other.coefficient)/(self.center-other.center).abs().powf(2.0));
        matrix_element *= (PI/(self.coefficient+other.coefficient)).powf(0.75);
        matrix_element *= -self.coefficient*other.coefficient/(self.coefficient+other.coefficient).exp() * (self.center-other.center).abs().powf(2.0);
        matrix_element

    }

    pub fn nuclear_attraction_integral(&self, other: &Gaussian, atom: &Atom) -> f32 {

        //Calculates the un-normalised nuclear attraction integrals
        let center_p = (self.coefficient*self.center + other.coefficient*other.center)/(self.coefficient + other.coefficient);

        let n = (2.0*self.coefficient/PI).powf(0.75) * (2.0*other.coefficient/PI).powf(0.75);
        let mut matrix_element = n*-2.0*PI/(self.coefficient+other.coefficient)*(atom.atomic_number as f32);
        matrix_element *= -self.coefficient*other.coefficient/(self.coefficient+other.coefficient).exp()*(self.center-other.center).abs().powf(2.0);

        if let Some(i) = Gaussian::f0((self.coefficient+other.coefficient)*(center_p-atom.position).abs().powf(2.0)) {
            matrix_element *= i;
        }
        matrix_element
    }
    
    pub fn two_electron_integral(&self, a: &Gaussian, b: &Gaussian, c: &Gaussian) -> f32 {

        /*
        Calculate two electron integrals
        self.coefficient, other.coefficient, c.coefficient, b.coefficient are the coefficients of the gaussian orbitals
        r_p is the distance between a and b, r_q is the distance between c and d
        */

        let center_p = (self.coefficient*self.center + a.coefficient*a.center)/(self.coefficient + a.coefficient);
        let center_q = (b.coefficient*b.center + c.coefficient*c.center)/(b.coefficient + c.coefficient);

        let n = (2.0*self.coefficient/PI).powf(0.75) * (2.0*a.coefficient/PI).powf(0.75)*(2.0*b.coefficient/PI).powf(0.75) * (2.0*c.coefficient/PI).powf(0.75);

        let mut matrix_element = n*2.0*PI.powf(2.5);
        matrix_element /= (self.coefficient+a.coefficient)*(b.coefficient+c.coefficient)*(self.coefficient+a.coefficient+b.coefficient+c.coefficient).sqrt();
        matrix_element *= -self.coefficient*a.coefficient/(self.coefficient+a.coefficient).exp()*(self.center-a.center).abs().powf(2.0) - b.coefficient*c.coefficient/(b.coefficient+c.coefficient)*(b.center-c.center).powf(2.0);

        if let Some(i) = Gaussian::f0((self.coefficient+a.coefficient)*(b.coefficient+c.coefficient)/(self.coefficient+a.coefficient+b.coefficient+c.coefficient)*(center_p-center_q).abs().powf(2.0)) {
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
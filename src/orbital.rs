use itertools::iproduct;
use serde::{Serialize, Deserialize};

use crate::{molecule::Atom, gaussian::Gaussian};

//Slater Type Orbital fit with N primative gausians (STO-NG) type basis
#[derive(Debug, Serialize, Deserialize)]
pub struct Orbital {
    //contraction coefficients
    pub exponents: Vec<f32>,

    //primitive Gaussians
    pub gaussians: Vec<Gaussian>
}

impl Orbital {

    pub fn two_center_contraction_with_atom<F>(&self, other: &Orbital, sto_ng: usize, atom: &Atom, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Atom) -> f32 {

        iproduct!(0..sto_ng, 0..sto_ng).fold(0.0, |acc, (i, j)| {
            acc + self.exponents[i]*other.exponents[j]*integral(&self.gaussians[i], &self.gaussians[j], atom)
        })
    }

    pub fn two_center_contraction<F>(&self, other: &Orbital, sto_ng: usize, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian) -> f32 {

        iproduct!(0..sto_ng, 0..sto_ng).fold(0.0, |acc, (i, j)| {
            acc + self.exponents[i]*other.exponents[j]*integral(&self.gaussians[i], &self.gaussians[j])
        })
    }

    pub fn four_center_contraction<F>(&self, a: &Orbital, b: &Orbital, c: &Orbital, sto_ng: usize, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Gaussian, &Gaussian) -> f32 {

        iproduct!(0..sto_ng, 0..sto_ng, 0..sto_ng, 0..sto_ng).fold(0.0, |acc, (i, j, k, l)| {
            acc + self.exponents[i]*a.exponents[j]*b.exponents[k]*c.exponents[l]
            * integral(&self.gaussians[i], &a.gaussians[j], &b.gaussians[k], &c.gaussians[l])
        })

    }
}
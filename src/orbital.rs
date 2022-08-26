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
        let mut total: f32 = Default::default();

        for i in 0..sto_ng {
            for j in 0..sto_ng {
                let gaussian_a = self.gaussians.get(i).unwrap();
                let gaussian_b = other.gaussians.get(j).unwrap();

                let exp_a = self.exponents.get(i).unwrap();
                let exp_b = other.exponents.get(j).unwrap();

                total += exp_a*exp_b*integral(gaussian_a, gaussian_b, atom)
            }
        }
        total
    }

    pub fn two_center_contraction<F>(&self, other: &Orbital, sto_ng: usize, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian) -> f32 {
        let mut total: f32 = Default::default();

        for i in 0..sto_ng {
            for j in 0..sto_ng {
                let gaussian_a = self.gaussians.get(i).unwrap();
                let gaussian_b = other.gaussians.get(j).unwrap();

                let exp_a = self.exponents.get(i).unwrap();
                let exp_b = other.exponents.get(j).unwrap();

                total += exp_a*exp_b*integral(gaussian_a, gaussian_b);
            }
        }
        total
    }

    pub fn four_center_contraction<F>(&self, a: &Orbital, b: &Orbital, c: &Orbital, sto_ng: usize, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Gaussian, &Gaussian) -> f32 {
        let mut total: f32 = Default::default();
        
        for i in 0..sto_ng {
            for j in 0..sto_ng {
                for k in 0..sto_ng {
                    for l in 0..sto_ng {

                        let gaussian_a = self.gaussians.get(i).unwrap();
                        let gaussian_b = a.gaussians.get(j).unwrap();
                        let gaussian_c = b.gaussians.get(k).unwrap();
                        let gaussian_d = c.gaussians.get(l).unwrap();

                        let exp_a = self.exponents.get(i).unwrap();
                        let exp_b = a.exponents.get(j).unwrap();
                        let exp_c = b.exponents.get(k).unwrap();
                        let exp_d = c.exponents.get(l).unwrap();

                        total += exp_a*exp_b*exp_c*exp_d*integral(gaussian_a, gaussian_b, gaussian_c, gaussian_d);
                    }
                }
            }
        }
        total
    }
}
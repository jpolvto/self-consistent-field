use serde::{Serialize, Deserialize};

use crate::{molecule::Atom, gaussian::Gaussian};

//Slater Type Orbital fit with N primative gausians (STO-NG) type basis
#[derive(Debug, Serialize, Deserialize)]
pub struct Orbital {
    //number n in stong
    pub n: usize,

    //contraction coefficients
    pub exponents: Vec<f32>,

    //primitive Gaussians
    pub gaussians: Vec<Gaussian>
}

impl Orbital {

    pub fn two_center_contraction_with_atom<F>(a: &Orbital, b: &Orbital, sto_ng: usize, atom: &Atom, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Atom) -> f32 {
        let mut total: f32 = Default::default();

        for i in 0..sto_ng {
            for j in 0..sto_ng {
                let gaussian_a = a.gaussians.get(i).unwrap();
                let gaussian_b = b.gaussians.get(j).unwrap();

                let exp_a = a.exponents.get(i).unwrap();
                let exp_b = b.exponents.get(j).unwrap();

                total += exp_a*exp_b*integral(gaussian_a, gaussian_b, atom)
            }
        }
        total
    }

    pub fn two_center_contraction<F>(a: &Orbital, b: &Orbital, sto_ng: usize, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian) -> f32 {
        let mut total: f32 = Default::default();

        for i in 0..sto_ng {
            for j in 0..sto_ng {
                let gaussian_a = a.gaussians.get(i).unwrap();
                let gaussian_b = b.gaussians.get(j).unwrap();

                let exp_a = a.exponents.get(i).unwrap();
                let exp_b = b.exponents.get(j).unwrap();

                total += exp_a*exp_b*integral(gaussian_a, gaussian_b);
            }
        }
        total
    }

    pub fn four_center_contraction<F>(a: &Orbital, b: &Orbital, c: &Orbital, d: &Orbital, sto_ng: usize, integral: F) -> f32
    where F: Fn(&Gaussian, &Gaussian, &Gaussian, &Gaussian) -> f32 {
        let mut total: f32 = Default::default();
        
        for i in 0..sto_ng {
            for j in 0..sto_ng {
                for k in 0..sto_ng {
                    for l in 0..sto_ng {

                        let gaussian_a = a.gaussians.get(i).unwrap();
                        let gaussian_b = b.gaussians.get(j).unwrap();
                        let gaussian_c = c.gaussians.get(k).unwrap();
                        let gaussian_d = d.gaussians.get(l).unwrap();

                        let exp_a = a.exponents.get(i).unwrap();
                        let exp_b = b.exponents.get(j).unwrap();
                        let exp_c = c.exponents.get(k).unwrap();
                        let exp_d = d.exponents.get(l).unwrap();

                        total += exp_a*exp_b*exp_c*exp_d*integral(gaussian_a, gaussian_b, gaussian_c, gaussian_d);
                    }
                }
            }
        }
        total
    }
}
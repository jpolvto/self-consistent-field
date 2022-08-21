use nalgebra::{Vector3, DMatrix};
use ndarray::Array4;
use serde::{Serialize, Deserialize};

use crate::{basis_set::BasisSet, gaussian::Gaussian};

#[derive(Serialize, Deserialize, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    // every atom must have an atomic number to identify the atom
    pub atomic_number: i32,

    //the position is specified as a 3D vector
    pub position: Vector3<f32>,
}

impl Molecule {
    pub fn nuclear_repulsion_energy(&self) -> f32 {
        let mut energy: f32 = 0.0;

        for atom_a in 0..self.atoms.len() {
            for atom_b in (atom_a + 1)..self.atoms.len() {
                energy += (self.atoms[atom_a].atomic_number) as f32 / (self.atoms[atom_a].position - self.atoms[atom_b].position).norm();
            }
        }
        energy
    }

    pub fn create_gaussians (&self, basis_set: &BasisSet, n: usize) -> Vec<Gaussian> {
        let mut gaussians: Vec<Gaussian> = Vec::new();

        for atom in &self.atoms {
            for shell in &basis_set.elements.get(&atom.atomic_number).unwrap().electron_shells {
                for orbital_coefficients in &shell.coefficients {
                    for i in 0..=(n-1) {
                        gaussians.push(Gaussian{ center: atom.position, coefficient: orbital_coefficients.get(i).unwrap().clone(), exponent: shell.exponents.get(i).unwrap().clone() })
                    }
                }
            }
        }
        gaussians
    }

    pub fn kinetic_matrix(gaussians: &Vec<Gaussian>, size: usize) -> DMatrix<f32> {
        let mut kinetic = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                let a  = gaussians.get(i).unwrap();
                let b = gaussians.get(j).unwrap();
                let (p, r_ab, k_ab) = Gaussian::gaussian_product(a, b);
                [kinetic[(i, j)], kinetic[(j, i)]] = [Gaussian::kinetic_energy_integral(a, b, &p, r_ab, k_ab); 2];
            }
        }
        kinetic
    }

    pub fn overlap_matrix(gaussians: &Vec<Gaussian>, size: usize) -> DMatrix<f32> {
        let mut overlap: DMatrix<f32> = DMatrix::identity(size, size);

        for i in 0..size {
            for j in (i+1)..size {
                let a  = gaussians.get(i).unwrap();
                let b = gaussians.get(j).unwrap();
                let (p, _r_ab, k_ab) = Gaussian::gaussian_product(a, b);
                [overlap[(i, j)], overlap[(j, i)]] = [Gaussian::overlap_integral(&p, k_ab); 2];
            }
        }
        overlap
    }

    pub fn nuclear_attraction_matrix(gaussians: &Vec<Gaussian>, size: usize, atoms: &Vec<Atom>) -> DMatrix<f32> {
        let mut nuclear_attraction = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                let a  = gaussians.get(i).unwrap();
                let b = gaussians.get(j).unwrap();
                let (p, _r_ab, k_ab) = Gaussian::gaussian_product(a, b);

                for atom in atoms {
                    let element = Gaussian::nuclear_attraction_integral(&atom, &p, k_ab);
                    nuclear_attraction[(i,j)] += element;
                    if i != j {
                        nuclear_attraction[(j,i)] += element;
                    }
                }
            }
        }
        nuclear_attraction    
    }

    pub fn two_electron_matrix(gaussians: &Vec<Gaussian>, size: usize) -> Array4<f32> {
        let mut two_electron = Array4::<f32>::zeros((size, size, size, size));

        for i in 0..size {
            for j in 0..size {
                for k in 0..size {
                    for l in 0..size {

                        let a  = gaussians.get(i).unwrap();
                        let b = gaussians.get(j).unwrap();
                        let c  = gaussians.get(k).unwrap();
                        let d = gaussians.get(l).unwrap();

                        let (p, _r_ab, k_ab) = Gaussian::gaussian_product(a, b);
                        let (q, _r_cd, k_cd) = Gaussian::gaussian_product(c, d);
                        [two_electron[[i, j, k, l]], two_electron[[i, k, l, j]]] = [Gaussian::two_electron_integral(&p, &q, k_ab, k_cd); 2]
                    }
                }
            }
        }
        two_electron
    }

    pub fn initial_values(overlap: DMatrix<f32>, kinetic: DMatrix<f32>, potential: DMatrix<f32>) -> (DMatrix<f32>, DMatrix<f32>) {

        let h_core = &kinetic*&potential;
        let overlap_eigen = overlap.symmetric_eigen();
        let mapped_eigenvalues = overlap_eigen.eigenvalues.map(|x: f32| x.powf(-0.5));
        let x = &overlap_eigen.eigenvectors * DMatrix::from_diagonal(&mapped_eigenvalues)*&overlap_eigen.eigenvectors.adjoint();

        (h_core, x)
    }

    pub fn hartree_fock(size: usize, h_core: DMatrix<f32>, nuclear_repulsion_energy: f32, x: &DMatrix<f32>, two_electron: Array4<f32>) -> (f32, f32) {

        let mut p = DMatrix::<f32>::zeros(size, size);
        let mut total_energy: f32 = Default::default();
        let mut old_energy:f32 = Default::default();
        let mut electronic_energy: f32 = Default::default();
        let treshold: f32 = 100.0;

        /*
         the term two_electron[[i, j, k, l]] is the coulomb coefficient
         the term two_electron[(i, k, l, j)] is the exchange coefficient
        */

        while treshold > 1e-03 {
            let mut g = DMatrix::<f32>::zeros(size, size);
            for i in 0..size {
                for j in 0..size {
                    for k in 0..size {
                        for l in k..size {
                            g[(i,j)] += p[(k,l)]*(two_electron[[i, j, k, l]]-0.5*two_electron[(i, k, l, j)]);
                        }
                    }
                }
            }

            let fock = &h_core + g;
            let f_prime = x.adjoint()*fock*x;
            let c = x*f_prime.symmetric_eigenvalues();

            for i in 0..size {
                for j in 0..size {
                    p[(i,j)] = 2.0*c[(i,0)]*c[(j,0)]
                }
            }

            electronic_energy += 100.0;

            total_energy = electronic_energy + nuclear_repulsion_energy;

            if total_energy - old_energy < 1e-6 {
                break
            }
    
            old_energy = total_energy;
        }


        (total_energy, electronic_energy)
    }
}

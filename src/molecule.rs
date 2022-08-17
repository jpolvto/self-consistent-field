use std::default;

use nalgebra::{Vector3, DMatrix, SquareMatrix, ComplexField};
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

    //open shell calculations
    number_of_electrons: i32,

    //open shell calculations
    principal_quantum_number: i32
}

impl Molecule {
    pub fn get_nuclear_repulsion_energy(&self) -> f32 {
        let mut energy: f32 = 0.0;

        for atom_a in 0..self.atoms.len() {
            for atom_b in (atom_a + 1)..self.atoms.len() {
                energy += (self.atoms[atom_a].atomic_number) as f32 / (self.atoms[atom_a].position - self.atoms[atom_b].position).norm();
            } 
        }
        energy;
    }

    pub fn get_gaussians (&self, basis_set: &BasisSet, n: usize) -> Vec<Gaussian> {
        let mut gaussians: Vec<Gaussian> = Vec::new();

        for atom in &self.atoms {
            for shell in &basis_set.elements.get(&atom.atomic_number).unwrap().electron_shells {
                for orbital_coefficients in &shell.coefficients {
                    for i in 1..n {
                        gaussians.push(Gaussian{ center: atom.position, coefficient: orbital_coefficients.get(n).unwrap().clone(), exponent: shell.exponents.get(n).unwrap().clone() })
                    }
                }
            }
        }
        gaussians
    }

    pub fn get_initial_values(size: usize, gaussians: &Vec<Gaussian>, molecule: Molecule) -> (DMatrix<f32>, DMatrix<f32>, DMatrix<f32>) {

        let mut kinetic = DMatrix::<f32>::zeros(size, size);
        let mut overlap = DMatrix::<f32>::zeros(size, size);
        let mut potential = DMatrix::<f32>::zeros(size, size);
        let mut multi_electron = DMatrix::<f32>::zeros(size, size, size, size);

        for i in 0..size {
            for j in i+1..size {

                let a  = gaussians.get(i).unwrap();
                let b = gaussians.get(j).unwrap();
                let (p, r_ab, k_ab) = Gaussian::gaussian_product(a, b);

                [kinetic[(i, j)], kinetic[(j, i)]] = [Gaussian::get_kinetic_integral(a, b, &p, r_ab, k_ab); 2];
                [overlap[(i, j)], overlap[(j, i)]] = [Gaussian::get_overlap_integral(&p, k_ab); 2];

                for atom in molecule.atoms {
                    [potential[(i,j)], potential[(j,i)]] = [Gaussian::get_potential_integral(&atom, &p, k_ab); 2];
                }

                for k in j+1..size {
                    for l in k+1..size {
                        let c  = gaussians.get(k).unwrap();
                        let d = gaussians.get(l).unwrap();

                        let (q, r_cd, k_cd) = Gaussian::gaussian_product(c, d);
                        [multi_electron[(i, j, k, l)], multi_electron[(i, k, l, j)]] = [Gaussian::get_multi_electron_integral(&p, &q, k_ab, k_cd); 2]
                    }
                }
            }
        }

        let h_core = kinetic+potential;
        let (s, u) = overlap.eigen_qr().unwrap();
        let x = u.dot((-1/2).exp())*u.adjoint();
        //X = U*diagm(s.^(-1/2))*U'   

        return (h_core, x, multi_electron)

    }

    pub fn hartree_fock(size: usize, gaussians: &Vec<Gaussian>, h_core: DMatrix<f32>, molecule: Molecule, x: DMatrix<f32>, multi_electron: DMatrix<f32>) -> (f32, f32) {
        let mut p_matrix = DMatrix::<f32>::zeros(size, size);

        let mut total_energy: f32 = Default::default();
        let mut old_energy:f32 = Default::default();
        let mut electronic_energy: f32 = Default::default();

        for scf_iteration in 0..100 {
            let mut g = DMatrix::<f32>::zeros(size, size);
            for i in 0..size {
                for j in i+1..size {
                    for k in j+1..size {
                        for l in k+1..size {
                            let a  = gaussians.get(i).unwrap();
                            let b = gaussians.get(j).unwrap();
                            let c  = gaussians.get(k).unwrap();
                            let d = gaussians.get(l).unwrap();

                            let (p, r_ab, k_ab) = Gaussian::gaussian_product(a, b);
                            let (q, r_cd, k_cd) = Gaussian::gaussian_product(c, d);

                            let coulomb  = multi_electron[(i, j, k, l)];
                            let exchange = multi_electron[(i, k, l, j)];

                            g[(i,j)] += p_matrix[(k,l)]*(coulomb-0.5*exchange);
                        }
                    }
                }
            }
            let f = h_core + g;
            
            let electronic_energy: f32 = Default::default();

            for i in 0..size {
                for j in 0..size {
                    electronic_energy += p[(j,i)]*(h_core[(i,j)]+f[(i,j)])
                }
            }

            total_energy = (electronic_energy*0.5) + Molecule::get_nuclear_repulsion_energy(&molecule);

            if scf_iteration > 1 && (old_energy - total_energy).abs() < 1e-6 {
                break
            }

            let f_prime = x.adjoint()*f*x;
            let (epsilon, c_prime) = f_prime.eigen_qr().unwrap();
            let c = (x*c_prime).real();

            let mut p = DMatrix::<f32>::zeros(size, size);

            for i in 0..size {
                for j in 0..size {
                    p[(i,j)] = 2.0*c[(i,0)]*c[(j,0)]
                }
            }
    
            old_energy = total_energy
        }
        (total_energy, electronic_energy)
    }
}

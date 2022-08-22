use nalgebra::{Vector3, DMatrix};
use ndarray::Array4;
use serde::{Serialize, Deserialize};

use crate::{basis_set::BasisSet, gaussian::{Orbital, Gaussian}};

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
    pub fn nuclear_attraction_energy(&self) -> f32 {
        let mut energy: f32 = 0.0;

        for atom_a in 0..self.atoms.len() {
            for atom_b in (atom_a + 1)..self.atoms.len() {
                energy += (self.atoms[atom_a].atomic_number) as f32 / (self.atoms[atom_a].position - self.atoms[atom_b].position).norm();
            }
        }
        energy
    }

    pub fn create_orbitals (&self, basisset: &BasisSet) -> Vec<Orbital> {
        let mut result: Vec<Orbital> = Vec::new();
        
        for atom in &self.atoms {
            for shell in &basisset.elements.get(&atom.atomic_number).unwrap().electron_shells {
                for orbital_coefficients in &shell.coefficients {
                    result.push(Orbital::create_orbital(orbital_coefficients, atom.position, &shell.exponents));
                }
            }
        }

        result
    }

    pub fn kinetic_matrix(orbitals: &Vec<Orbital>, size: usize) -> DMatrix<f32> {
        let mut kinetic = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                let a  = orbitals.get(i).unwrap();
                let b = orbitals.get(j).unwrap();

                [kinetic[(i, j)], kinetic[(j, i)]] = [Orbital::two_center_contraction(a, b, &Atom { atomic_number: 0, position: Vector3::new(0.0, 0.0, 0.0) },Gaussian::kinetic_energy_integral); 2];
            }
        }
        kinetic
    }

    pub fn overlap_matrix(orbitals: &Vec<Orbital>, size: usize) -> DMatrix<f32> {
        let mut overlap: DMatrix<f32> = DMatrix::identity(size, size);

        for i in 0..size {
            for j in (i+1)..size {
                let a  = orbitals.get(i).unwrap();
                let b = orbitals.get(j).unwrap();

                [overlap[(i, j)], overlap[(j, i)]] = [Orbital::two_center_contraction(a, b, &Atom { atomic_number: 0, position: Vector3::new(0.0, 0.0, 0.0) },Gaussian::overlap_integral); 2];
            }
        }
        overlap
    }

    pub fn nuclear_attraction_matrix(orbitals: &Vec<Orbital>, size: usize, atoms: &Vec<Atom>) -> DMatrix<f32> {
        let mut nuclear_attraction = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                let a  = orbitals.get(i).unwrap();
                let b = orbitals.get(j).unwrap();

                for atom in atoms {
                    let element = Orbital::two_center_contraction(a, b, &Atom { atomic_number: 0, position: Vector3::new(0.0, 0.0, 0.0) },Gaussian::nuclear_attraction_integral);
                    nuclear_attraction[(i,j)] += element;
                    if i != j {
                        nuclear_attraction[(j,i)] += element;
                    }
                }
            }
        }
        nuclear_attraction    
    }

    pub fn two_electron_matrix(orbitals: &Vec<Orbital>, size: usize) -> Array4<f32> {
        let mut two_electron = Array4::<f32>::zeros((size, size, size, size));

        for i in 0..size {
            for j in 0..size {
                for k in 0..size {
                    for l in 0..size {

                        let a  = orbitals.get(i).unwrap();
                        let b = orbitals.get(j).unwrap();
                        let c  = orbitals.get(k).unwrap();
                        let d = orbitals.get(l).unwrap();
                    
                        [two_electron[[i, j, k, l]], two_electron[[i, k, l, j]]] = [Orbital::four_center_contraction(a, b, c, d, Gaussian::two_electron_integral); 2]
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

    pub fn hartree_fock(size: usize, h_core: DMatrix<f32>, nuclear_attraction_energy: f32, x: &DMatrix<f32>, two_electron: Array4<f32>) -> (f32, f32, i32) {

        let mut old_energy:f32 = Default::default();
        let mut electronic_energy: f32 = Default::default();
        let scf_max = 1000;
        let mut iterations = Default::default();

        /*
         the term two_electron[[i, j, k, l]] is the coulomb coefficient
         the term two_electron[(i, k, l, j)] is the exchange coefficient
        */

        for scf in 0..scf_max {
            let mut total_energy: f32 = Default::default();
            let mut p = DMatrix::<f32>::zeros(size, size);
            let mut g = DMatrix::<f32>::zeros(size, size);
            
            for i in 0..size {
                for j in 0..size {
                    for k in 0..size {
                        for l in 0..size {
                            g[(i,j)] += p[(k,l)]*(two_electron[[i, j, k, l]]-0.5*two_electron[(i, k, l, j)]);
                        }
                    }
                }
            }

            let fock = &h_core + g;

            for i in 0..size {
                for j in 0..size {
                    electronic_energy += p[(i,0)]*h_core[(i,j)]+fock[(i,j)]
                }
            }

            total_energy = (electronic_energy*0.5) + nuclear_attraction_energy;

            if (old_energy - total_energy).abs() < 1e-6 {
                iterations = scf;
                break
            }

            let f_prime = x.adjoint()*fock*x;
            let c = x*f_prime.symmetric_eigenvalues();

            for i in 0..size {
                for j in 0..size {
                    p[(i,j)] = 2.0*c[(i,0)]*c[(j,0)]
                }
            }

            old_energy = total_energy;
        }
        (old_energy, electronic_energy, iterations)
    }
}

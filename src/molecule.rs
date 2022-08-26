use nalgebra::{Vector3, DMatrix};
use ndarray::Array4;
use serde::{Serialize, Deserialize};
use crate::gaussian::Gaussian;
use crate::orbital::Orbital;

#[derive(Serialize, Deserialize, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub orbitals: Vec<Orbital>,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub sto_ng: usize
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Atom {
    // every atom must have an atomic number to identify the atom
    pub atomic_number: usize,

    //the position is specified as a 3D vector
    pub position: Vector3<f32>,
}

impl Molecule {
    pub fn nuclear_attraction_energy(&self) -> f32 {
        let mut energy: f32 = 0.0;

        for atom_a in 0..self.atoms.len() {
            for atom_b in (atom_a + 1)..self.atoms.len() {
                energy += (self.atoms[atom_a].atomic_number as f32) / (self.atoms[atom_a].position - self.atoms[atom_b].position).norm();
            }
        }
        energy
    }

    pub fn create_orbitals(&mut self) {

        let mut orbitals: Vec<Orbital> = Default::default();
        
        for atom in &self.atoms {
            let center = atom.position;

            orbitals.push(Orbital {
                exponents: Vec::from([0.444635, 0.535328, 0.154329]),
                gaussians: Vec::from([Gaussian { center, coefficient: 0.444635 }, Gaussian { center, coefficient: 0.535328 }, Gaussian { center, coefficient: 0.154329 }])
            })
        }

        self.orbitals = orbitals;
        self.sto_ng = 3;
    }

    pub fn kinetic_matrix(&self) -> DMatrix<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let sto_ng = self.sto_ng;

        let mut kinetic = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                let a  = orbitals.get(i).unwrap();
                let b = orbitals.get(j).unwrap();

                [kinetic[(i, j)], kinetic[(j, i)]] = [a.two_center_contraction(b, sto_ng, Gaussian::kinetic_energy_integral); 2];
                
            }
        }
        kinetic
    }

    pub fn overlap_matrix(&self) -> DMatrix<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let sto_ng = self.sto_ng;
        let mut overlap: DMatrix<f32> = DMatrix::identity(size, size);

        for i in 0..size {
            for j in (i+1)..size {
                let a  = orbitals.get(i).unwrap();
                let b = orbitals.get(j).unwrap();

                [overlap[(i, j)], overlap[(j, i)]] = [a.two_center_contraction(b, sto_ng, Gaussian::overlap_integral); 2];
            }
        }
        overlap
    }

    pub fn nuclear_attraction_matrix(&self) -> DMatrix<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let sto_ng = self.sto_ng;
        let mut nuclear_attraction = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                let a  = orbitals.get(i).unwrap();
                let b = orbitals.get(j).unwrap();

                for atom in &self.atoms {
                    let element = a.two_center_contraction_with_atom(b, sto_ng, atom, Gaussian::nuclear_attraction_integral);
                    nuclear_attraction[(i,j)] += element;
                    if i != j {
                        nuclear_attraction[(j,i)] += element;
                    }
                }
            }
        }
        nuclear_attraction    
    }

    pub fn two_electron_matrix(&self) -> Array4<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let sto_ng = self.sto_ng;
        let mut two_electron = Array4::<f32>::zeros((size, size, size, size));

        for i in 0..size {
            for j in 0..size {
                for k in 0..size {
                    for l in 0..size {

                        let a  = orbitals.get(i).unwrap();
                        let b = orbitals.get(j).unwrap();
                        let c  = orbitals.get(k).unwrap();
                        let d = orbitals.get(l).unwrap();
                    
                        [two_electron[[i, j, k, l]], two_electron[[i, k, l, j]]] = [a.four_center_contraction(b, c, d, sto_ng, Gaussian::two_electron_integral); 2]
                    }
                }
            }
        }
        two_electron
    }

    pub fn hartree_fock(&self) {

        let two_electron = self.two_electron_matrix();
        let kinetic = self.kinetic_matrix();
        let overlap = self.overlap_matrix();
        let nuclear_attraction_matrix = self.nuclear_attraction_matrix();
        let nuclear_attraction_energy = self.nuclear_attraction_energy();
        let size = self.orbitals.len();

        let mut total_energy:f32 = Default::default();
        let mut electronic_energy: f32 = Default::default();
        let mut p = DMatrix::<f32>::zeros(size, size);
        let mut p_previous = DMatrix::<f32>::zeros(size, size);
        let mut p_list: Vec<DMatrix<f32>> = Default::default();

        let h_core = &kinetic+&nuclear_attraction_matrix;
        let overlap_eigen = overlap.symmetric_eigen();
        let mapped_eigenvalues = overlap_eigen.eigenvalues.map(|x: f32| x.powf(-0.5));
        let x = &overlap_eigen.eigenvectors * DMatrix::from_diagonal(&mapped_eigenvalues)*&overlap_eigen.eigenvectors.adjoint();

        /*
         the term two_electron[[i, j, k, l]] is the coulomb coefficient
         the term two_electron[(i, k, l, j)] is the exchange coefficient
        */

        for scf_iter in 0..100 {

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

            let f = &h_core + g;

            for i in 0..size {
                for j in 0..size {
                    electronic_energy += p[(i,0)]*h_core[(i,j)]+f[(i,j)]
                }
            }

            let f_prime = &x.adjoint()*f*&x;
            let f_prime_eigen = f_prime.symmetric_eigen();
            let c_prime = &f_prime_eigen.eigenvectors;
            let c = &x*c_prime;

            for i in 0..size {
                for j in 0..size {
                    p[(i,j)] = 2.0*c[(i,0)]*c[(j,0)]
                }
            }

        }
        println!("Self consistent field finished with total energy: {} ", total_energy);
    }
}
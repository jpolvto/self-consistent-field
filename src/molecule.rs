use itertools::iproduct;
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
        
        self.orbitals = self.atoms.iter().map(|a| {
            Orbital {
                exponents: Vec::from([0.444635, 0.535328, 0.154329]),
                gaussians: Vec::from([Gaussian { center: a.position, coefficient: 0.444635 }, Gaussian { center: a.position, coefficient: 0.535328 }, Gaussian { center: a.position, coefficient: 0.154329 }])
            }
        }).collect();
        
        self.sto_ng = 3;
    }

    pub fn kinetic_matrix(&self) -> DMatrix<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let mut kinetic = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                [kinetic[(i, j)], kinetic[(j, i)]] = [orbitals[i].two_center_contraction(&orbitals[j], self.sto_ng, Gaussian::kinetic_energy_integral); 2]; 
            }
        }
        kinetic
    }

    pub fn overlap_matrix(&self) -> DMatrix<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let mut overlap: DMatrix<f32> = DMatrix::identity(size, size);

        for i in 0..size {
            for j in (i+1)..size {
                [overlap[(i, j)], overlap[(j, i)]] = [orbitals[i].two_center_contraction(&orbitals[j], self.sto_ng, Gaussian::overlap_integral); 2];
            }
        }
        overlap
    }

    pub fn nuclear_attraction_matrix(&self) -> DMatrix<f32> {
        let orbitals = &self.orbitals;
        let size = orbitals.len();
        let mut nuclear_attraction = DMatrix::<f32>::zeros(size, size);

        for i in 0..size {
            for j in i..size {
                for atom in &self.atoms {
                    let element = orbitals[i].two_center_contraction_with_atom(&orbitals[j], self.sto_ng, atom, Gaussian::nuclear_attraction_integral);
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
        let mut two_electron = Array4::<f32>::zeros((size, size, size, size));

        iproduct!(0..size, 0..size, 0..size, 0..size).for_each(|(i, j, k, l)| {        
            [two_electron[[i, j, k, l]], two_electron[[i, k, l, j]]] = [orbitals[i].four_center_contraction(&orbitals[j],& orbitals[k], &orbitals[l], self.sto_ng, Gaussian::two_electron_integral); 2]            
        });
        two_electron
    }
}
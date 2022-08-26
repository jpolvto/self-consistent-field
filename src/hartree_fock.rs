use std::{env, fs};
use nalgebra::DMatrix;

use crate::molecule::Molecule;

pub fn hartree_fock() {

    dotenv::dotenv().ok();
    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let mut molecule: Molecule = serde_json::from_str(&input_contents).unwrap();
    molecule.create_orbitals();

    let two_electron = molecule.two_electron_matrix();
    let kinetic = molecule.kinetic_matrix();
    let overlap = molecule.overlap_matrix();
    let nuclear_attraction_matrix = molecule.nuclear_attraction_matrix();
    let nuclear_attraction_energy = molecule.nuclear_attraction_energy();
    let size = molecule.orbitals.len();

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
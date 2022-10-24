use std::{env, fs};
use itertools::iproduct;
use nalgebra::{DMatrix};

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

    println!("two_electron: {two_electron}, kinetic: {kinetic}, overlap: {overlap}, nuclear_attraction_matrix: {nuclear_attraction_matrix}, nuclear_attraction_energy: {nuclear_attraction_energy}");

    let mut old_energy: f32 = Default::default();
    let mut total_energy: f32 = Default::default();
    let mut count: usize = Default::default();

    // calculating h_core using the kinetic and nuclear_attraction_matrix

    let h_core = &kinetic+&nuclear_attraction_matrix;
    let overlap_eigen = overlap.symmetric_eigen();
    let mapped_eigenvalues = overlap_eigen.eigenvalues.map(|x: f32| x.powf(-0.5));
    let x = &overlap_eigen.eigenvectors * DMatrix::from_diagonal(&mapped_eigenvalues)*&overlap_eigen.eigenvectors.adjoint();

    let mut p = DMatrix::<f32>::zeros(size, size);

    loop {

        count += 1;
        let mut g = DMatrix::<f32>::zeros(size, size);

        /*
        the term two_electron[[i, j, k, l]] is the coulomb coefficient
        the term two_electron[(i, k, l, j)] is the exchange coefficient
        */

        iproduct!(0..size, 0..size, 0..size, 0..size).for_each(|(i, j, k, l)| {
            g[(i,j)] += p[(k,l)]*(two_electron[[i, j, k, l]]-0.5*two_electron[(i, k, l, j)]);            
        });
        
    }
    println!("Self consistent field finished with total energy: {} ", total_energy);
}

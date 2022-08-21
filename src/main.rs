mod molecule;
mod basis_set;
mod gaussian;

use std::env;
use std::fs;

use basis_set::BasisSet;
use gaussian::Gaussian;
use molecule::Molecule;
fn main() {
    dotenv::dotenv().ok();

    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let molecule: Molecule = serde_json::from_str(&input_contents).unwrap();

    let basis_set_file_name = env::var("SELF_CONSISTENT_FIELD_BASIS_SET_FILE_NAME").unwrap();
    let basis_set_contents = fs::read_to_string(basis_set_file_name).unwrap();
    let basis_set: BasisSet = serde_json::from_str(&basis_set_contents).unwrap();

    let gaussians: Vec<Gaussian> = molecule.create_gaussians(&basis_set, 3);

    //println!("gaussians: {}", serde_json::to_string(&gaussians).unwrap());

    let size = gaussians.len() -1;

    //println!("size: {}", size);

    let nuclear_repulsion_energy = molecule.nuclear_repulsion_energy();

    //println!("nuclear_repulsion_energy: {}", nuclear_repulsion_energy);

    let atoms = molecule.atoms;

    //println!("atoms: {}", serde_json::to_string(&atoms).unwrap());

    let two_electron = Molecule::two_electron_matrix(&gaussians, size);

    //println!("two_electron: {}", two_electron);

    let kinetic = Molecule::kinetic_matrix(&gaussians, size);

    //println!("kinetic: {}", kinetic);

    let overlap = Molecule::overlap_matrix(&gaussians, size);

    //println!("overlap: {}", overlap);

    let nuclear_attraction_matrix = Molecule::nuclear_attraction_matrix(&gaussians, size, &atoms);

    //println!("nuclear_attraction_matrix: {}", nuclear_attraction_matrix);

    let (h_core, x) = Molecule::initial_values(overlap, kinetic, nuclear_attraction_matrix);

    //println!("h_core: {}\nx: {}", h_core, x);

    let (total_energy, electronic_energy, iterations) = Molecule::hartree_fock(size, h_core, nuclear_repulsion_energy, &x, two_electron);

    println!("self consistent field finished in {} iterations\ntotal energy: {}\nelectronic energy: {} ", iterations, total_energy, electronic_energy);
}

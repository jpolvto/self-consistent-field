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

    let gaussians: Vec<Gaussian> = molecule.get_gaussians(&basis_set, 3);
    let size = gaussians.len();

    let (h_core, x, multi_electron) = Molecule::get_initial_values(size, &gaussians, molecule);
    let (total_energy, electronic_energy) = Molecule::hartree_fock(size, &gaussians, h_core, molecule, x, multi_electron);

}

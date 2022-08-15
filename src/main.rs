mod molecule;
mod basis_set;
mod gaussian;

use std::env;
use std::fs;

use molecule::Molecule;
fn main() {
    dotenv::dotenv().ok();

    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let molecule: molecule::Molecule = serde_json::from_str(&input_contents).unwrap();

    let basis_set_file_name = env::var("SELF_CONSISTENT_FIELD_BASIS_SET_FILE_NAME").unwrap();
    let basis_set_contents = fs::read_to_string(basis_set_file_name).unwrap();
    let basis_set: basis_set::BasisSet = serde_json::from_str(&basis_set_contents).unwrap();

    let wavefunction = molecule.get_slater_orbitals(&basis_set, 3);

    let nuclear_repulsion_energy = molecule.get_nuclear_repulsion_energy();
}

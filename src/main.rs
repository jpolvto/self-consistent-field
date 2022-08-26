mod molecule;
mod orbital;
mod gaussian;

use std::env;
use std::fs;

use molecule::Molecule;

fn main() {
    dotenv::dotenv().ok();
    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let mut molecule: Molecule = serde_json::from_str(&input_contents).unwrap();
    molecule.create_orbitals(3);
    molecule.hartree_fock();

}

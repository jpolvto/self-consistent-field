mod molecule;
mod basisset;

use std::env;
use std::fs;


fn main() {
    dotenv::dotenv().ok();

    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let input_deserialized: molecule::Molecule = serde_json::from_str(&input_contents).unwrap();
    let input_serialized = serde_json::to_string(&input_deserialized).unwrap();

    let basis_set_file_name = env::var("SELF_CONSISTENT_FIELD_BASIS_SET_FILE_NAME").unwrap();
    let basis_set_contents = fs::read_to_string(basis_set_file_name).unwrap();
    let basis_set_deserialized: basisset::BasisSet = serde_json::from_str(&basis_set_contents).unwrap();
    let basis_set_serialized = serde_json::to_string(&basis_set_deserialized).unwrap();
    

    println!("Molecule:\n{}", input_serialized);
    println!("Molecule:\n{}", basis_set_serialized);
}

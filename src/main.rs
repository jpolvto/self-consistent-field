mod molecule;

use std::env;
use std::fs;


fn main() {
    dotenv::dotenv().ok();

    let file_name = env::var("SELF_CONSISTENT_FIELD_FILE_NAME").unwrap();
    let contents = fs::read_to_string(file_name).unwrap();
    let deserialized: molecule::Molecule = serde_json::from_str(&contents).unwrap();

    let bla = serde_json::to_string(&deserialized).unwrap();

    println!("Molecule:\n{}", bla);
}

mod molecule;

use std::env;
use std::fs;

use nalgebra::DMatrix;


fn main() {
    dotenv::dotenv().ok();

    let file_name = env::var("SELF_CONSISTENT_FIELD_FILE_NAME").unwrap();
    let contents = fs::read_to_string(file_name).unwrap();
    let deserialized: DMatrix<f32> = serde_json::from_str(&contents).unwrap();
    let bla = deserialized.eigenvalues().unwrap();

    // Compare against reference
    println!("Computed eigenvalues:\n{}", bla);
}

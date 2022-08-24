mod molecule;
mod orbital;
mod gaussian;

use std::env;
use std::fs;

use molecule::Molecule;

use crate::orbital::Orbital;


fn main() {
    dotenv::dotenv().ok();

    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let molecule: Molecule = serde_json::from_str(&input_contents).unwrap();

    const STO_NG: usize = 3;

    let orbitals: Vec<Orbital> = molecule.create_orbitals();

    println!("orbitals: {}", serde_json::to_string(&orbitals).unwrap());

    let size = orbitals.len();

    println!("size: {}", size);

    let nuclear_attraction_energy = molecule.nuclear_attraction_energy();

    println!("nuclear_attraction_energy: {}", nuclear_attraction_energy);

    let atoms = molecule.atoms;

    println!("atoms: {}", serde_json::to_string(&atoms).unwrap());

    let two_electron = Molecule::two_electron_matrix(&orbitals, STO_NG, size);

    println!("two_electron: {}", two_electron);

    let kinetic = Molecule::kinetic_matrix(&orbitals, STO_NG, size);

    println!("kinetic: {}", kinetic);

    let overlap = Molecule::overlap_matrix(&orbitals, STO_NG, size);

    println!("overlap: {}", overlap);

    let nuclear_attraction_matrix = Molecule::nuclear_attraction_matrix(&orbitals, STO_NG, size, &atoms);

    println!("nuclear_attraction_matrix: {}", nuclear_attraction_matrix);

    let (total_energy, electronic_energy, iterations) = Molecule::hartree_fock(size, overlap, kinetic, nuclear_attraction_matrix, two_electron, nuclear_attraction_energy);

    println!("self consistent field finished in {} iterations\ntotal energy: {}\nelectronic energy: {} ", iterations, total_energy, electronic_energy);
}

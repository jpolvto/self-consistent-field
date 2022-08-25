mod molecule;
mod orbital;
mod gaussian;

use std::env;
use std::fs;

use molecule::Molecule;
use orbital::Orbital;

const STO_NG: usize = 3;

fn main() {
    dotenv::dotenv().ok();

    let input_file_name = env::var("SELF_CONSISTENT_FIELD_INPUT_FILE_NAME").unwrap();
    let input_contents = fs::read_to_string(input_file_name).unwrap();
    let molecule: Molecule = serde_json::from_str(&input_contents).unwrap();
    let orbitals: Vec<Orbital> = molecule.create_orbitals();

    let nuclear_attraction_energy = molecule.nuclear_attraction_energy();
    let two_electron = Molecule::two_electron_matrix(&orbitals, STO_NG);
    let kinetic = Molecule::kinetic_matrix(&orbitals, STO_NG);
    let overlap = Molecule::overlap_matrix(&orbitals, STO_NG);
    let nuclear_attraction_matrix = Molecule::nuclear_attraction_matrix(&orbitals, STO_NG, &molecule.atoms);

    //println!("nuclear_attraction_energy: {}", nuclear_attraction_energy);
    //println!("molecule: {}", serde_json::to_string(&molecule).unwrap());
    //println!("two_electron: {}", two_electron);
    //println!("kinetic: {}", kinetic);
    //println!("overlap: {}", overlap);
    //println!("orbitals: {}", serde_json::to_string(&orbitals).unwrap());
    //println!("nuclear_attraction_matrix: {}", nuclear_attraction_matrix);

    let (total_energy, electronic_energy, iterations) = Molecule::hartree_fock(orbitals.len(), overlap, kinetic, nuclear_attraction_matrix, two_electron, nuclear_attraction_energy);
    println!("self consistent field finished in {} iterations\ntotal energy: {}\nelectronic energy: {} ", iterations, total_energy, electronic_energy);
}

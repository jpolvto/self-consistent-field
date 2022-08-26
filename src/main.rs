mod molecule;
mod orbital;
mod gaussian;
mod hartree_fock;

use crate::hartree_fock::hartree_fock;

use std::env;
use std::fs;

use molecule::Molecule;

fn main() {
    hartree_fock();
}


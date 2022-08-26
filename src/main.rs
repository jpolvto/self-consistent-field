mod molecule;
mod orbital;
mod gaussian;
mod hartree_fock;

use crate::hartree_fock::hartree_fock;

fn main() {
    hartree_fock();
}


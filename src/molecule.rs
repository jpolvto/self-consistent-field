use nalgebra::Vector3;

struct Atom {
    id: i32,
    z: i32,
    position: Vector3<i32>,
    number_of_electrons: i32,
}

struct Molecule {
    atoms: Vec<AtomWithShells>,
    alpha_electrons: i32,
    beta_electrons: i32,
    count_number_of_contracted_gaussians: i32,
    count_number_of_gaussians: i32,
    nuclear_repulsion_energy: f32,
    number_of_electrons: i32,
    get_max_angular_momentum: i32,
}
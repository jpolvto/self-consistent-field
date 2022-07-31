use nalgebra::Vector3;

struct Atom {
    id: i32,
    z: i32,
    position: Vector3<i32>,
    number_of_electrons: i32,
}

struct Molecule {
    atoms: Vec<Atom>,
}

impl Molecule {
    fn get_max_angular_momentum(&self) {

    }

    fn get_number_of_contracted_gaussians(&self) {

    }

    fn get_number_of_gaussians(&self) {
        
    }

    fn get_alpha_electrons(&self) {
        
    }

    fn get_beta_electrons(&self) {
        
    }

    fn get_nuclear_repulsion_energy(&self) {
        
    }

    fn get_number_of_electrons(&self) {
        
    }

}
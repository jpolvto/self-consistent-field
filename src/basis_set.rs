use std::{collections::HashMap};
use serde::{Serialize, Deserialize};
use serde_with::{self, serde_as};
use serde_with::DisplayFromStr;

#[derive(Serialize, Deserialize, Debug)]
pub struct BasisSet {
    pub elements: HashMap<i32, Element>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Element {
    pub electron_shells: Vec<ElectronShell>
}

#[serde_as]
#[derive(Serialize, Deserialize, Debug)]
pub struct ElectronShell {
    #[serde_as(as = "Vec<DisplayFromStr>")]
    pub exponents: Vec<f32>,

    #[serde_as(as = "Vec<Vec<DisplayFromStr>>")]
    pub coefficients: Vec<Vec<f32>>
}

use std::collections::HashMap;
use std::hash::RandomState;
use std::rc::Rc;

use crate::domain::DomainContainer;
use crate::reference::Reference;
use crate::domain::Domain;


type RefRc = Rc<Reference>;
type StringHashMap<T> = HashMap<String, T, RandomState>;

#[derive(Debug)]
pub struct Alignment{
    matching_reference: Rc<Reference>,
    bitscore: f64
}

impl DomainContainer for Alignment {
    fn get_domains(&self) -> &Vec<Rc<Domain>> {
        &self.matching_reference.get_domains()
    }
    fn get_name(&self) -> &String {
        &self.matching_reference.get_name()
    }
}

impl Alignment {
    pub fn new(row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) -> Self {
        let matching_reference = ref_hashmap.get(row[1]).unwrap().clone();
        let bitscore = row[11].parse().unwrap();
        Self { matching_reference, bitscore}
    }
    pub fn get_bitscore(&self) -> &f64 {
        &self.bitscore
    }
}
use std::collections::HashMap;
use std::hash::RandomState;
use std::rc::Rc;

type StringHashMap<T> = HashMap<String, T, RandomState>;
type RcStringHashMap<T> = HashMap<Rc<String>, T, RandomState>;

#[derive(Debug)]
pub struct Domain {
    pub cdd_acc: String,
    pub super_acc: String
}
pub trait DomainContainer{
    fn get_domains(&self) -> &Vec<Rc<Domain>>;
    fn get_name(&self) -> &String;

    fn get_gene_name(&self) -> String {
        self.get_name().split("|").last().unwrap().to_string()
    }

    fn get_classification(&self) -> String {
        let fields: Vec<&str> = self.get_name().split("|").collect();
        fields[fields.len() - 2].to_string()
    }

    fn get_domain_identifiers(&self) -> (String, String, String) {
        if self.get_gene_name().to_uppercase() == "ARNA" && self.get_classification() == "peptide" {
            //print!("{}\n", self.get_name())
        }
        let gene_key = format!("{}|{}|{}", 
            self.get_domains().iter().map(|x| x.cdd_acc.clone()).collect::<Vec<_>>().join("$"), 
            self.get_gene_name().to_uppercase(), self.get_classification());
        let cdd_key = format!("{}|{}", 
            self.get_domains().iter().map(|x| x.cdd_acc.clone()).collect::<Vec<_>>().join("$"), 
            self.get_classification());
        let super_key = format!("{}|{}", 
            self.get_domains().iter().map(|x| x.super_acc.clone()).collect::<Vec<_>>().join("$"), 
            self.get_classification());
        (gene_key, cdd_key, super_key)
    }

    fn update_amr_hashmaps(&self, gene_map: &mut RcStringHashMap<usize>, 
            cdd_map: &mut RcStringHashMap<usize>, super_map: &mut RcStringHashMap<usize>,
            domain_name_refs: &mut StringHashMap<Rc<String>>){
        let (gene_key, cdd_key, super_key) = self.get_domain_identifiers();
        let gene_key_ref = {
            if domain_name_refs.contains_key(&gene_key) {
                domain_name_refs.get(&gene_key).unwrap().clone()
            }
            else {
                domain_name_refs.insert(gene_key.clone(), Rc::new(gene_key.clone()));
                domain_name_refs.get(&gene_key).unwrap().clone()
            }
        };
        if gene_map.contains_key(&gene_key_ref) {
            *gene_map.get_mut(&gene_key_ref).unwrap() += 1
            
        }
        else {
            gene_map.insert(gene_key_ref, 1);
        }
        let cdd_key_ref = {
            if domain_name_refs.contains_key(&cdd_key) {
                domain_name_refs.get(&cdd_key).unwrap().clone()
            }
            else {
                domain_name_refs.insert(cdd_key.clone(), Rc::new(cdd_key.clone()));
                domain_name_refs.get(&cdd_key).unwrap().clone()
            }
        };
        if cdd_map.contains_key(&cdd_key_ref) {
            *cdd_map.get_mut(&cdd_key_ref).unwrap() += 1;
        }
        else {
            cdd_map.insert(cdd_key_ref, 1);
        }
        let super_key_ref = {
            if domain_name_refs.contains_key(&super_key) {
                domain_name_refs.get(&super_key).unwrap().clone()
            }
            else {
                domain_name_refs.insert(super_key.clone(), Rc::new(super_key.clone()));
                domain_name_refs.get(&super_key).unwrap().clone()
            }
        };
        if super_map.contains_key(&super_key_ref) {
            *super_map.get_mut(&super_key_ref).unwrap() += 1;
        }
        else {
            super_map.insert(super_key_ref, 1);
        }
    }
}
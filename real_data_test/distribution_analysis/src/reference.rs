use crate::domain::DomainContainer;
use crate::domain::Domain;
use std::rc::Rc;

#[derive(Debug)]
pub struct Reference {
    name: String,
    domains: Vec<Rc<Domain>>
}

impl DomainContainer for Reference {
    fn get_domains(&self) -> &Vec<Rc<Domain>> {
        &self.domains
    }
    fn get_name(&self) -> &String {
        &self.name
    }
}

impl Reference {
    pub fn new(row: Vec<&str>) -> Self {
        let name: &str = row[0].split('>').last().unwrap();
        Self {
            name: name.to_string(),
            domains: vec![Rc::new(Domain { 
                cdd_acc: row[7].to_string(), 
                super_acc: if row[1] == "specific" {
                    row[10].to_string()
                } else {
                    row[7].to_string() 
                }
            })]
        }
    }
    pub fn new_raw(name: String, domains: Vec<Rc<Domain>>) -> Self{
        Self{name, domains}
    }
    pub fn add_to_domain(&mut self, row: Vec<&str>) {
        self.domains.push(Rc::new(Domain { 
            cdd_acc: row[7].to_string(), 
            super_acc: if row[1] == "specific" {
                row[10].to_string()
            } else {
                row[7].to_string() 
            }
        }))
    }
    
}
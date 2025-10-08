use std::collections::HashMap;
use std::hash::RandomState;
use std::rc::Rc;

use crate::domain::DomainContainer;
use crate::alignment::Alignment;
use crate::reference::Reference;

type RefRc = Rc<Reference>;
type StringHashMap<T> = HashMap<String, T, RandomState>;
type RcStringHashMap<T> = HashMap<Rc<String>, T, RandomState>;


#[derive(Debug)]
pub struct Query{
    alignments: Vec<Alignment>,
    top_diamond_alignment: usize,      //Should be the first alignment
    top_deeparg_hit: Option<usize>,    //Should be None until parsing DeepARG output
    is_deeparg_hit: bool
}

impl Query{
    pub fn new(row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) -> Self {
        Self { 
            alignments: vec![Alignment::new(row, ref_hashmap)],
            top_diamond_alignment: 0usize,
            top_deeparg_hit: None,
            is_deeparg_hit: false }
    }

    pub fn is_deeparg_hit(&self) -> &bool {
        &self.is_deeparg_hit
    }

    pub fn are_diamond_and_deeparg_the_same(&self) -> bool {
        if self.top_deeparg_hit.is_none() {
            panic!("You didn't check to see if query is a deeparg hit")
        }
        self.top_diamond_alignment == self.top_deeparg_hit.unwrap()
    }

    pub fn get_top_diamond_alignment_domain_identifiers(&self) -> (String, String, String) {
        self.alignments[self.top_diamond_alignment].get_domain_identifiers()
    }

    pub fn get_top_deeparg_hit_domain_identifiers(&self) -> Result<(String, String, String), &str> {
        if self.top_deeparg_hit.is_none() {
            Err("Query is not part of DeepARG output")
        }
        else {
            Ok(self.alignments[self.top_deeparg_hit.unwrap()].get_domain_identifiers())
        }
    }
    pub fn get_top_diamond_alignment(&self) -> &Alignment {
        &self.alignments[self.top_diamond_alignment]
    }
    pub fn get_top_deeparg_alignment(&self) -> &Alignment {
        &self.alignments[self.top_deeparg_hit.unwrap()]
    }

    pub fn add_alignment(&mut self, row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) {
        self.alignments.push(Alignment::new(row, ref_hashmap));
        if self.alignments.last().unwrap().get_bitscore() > self.alignments[self.top_diamond_alignment].get_bitscore() {
            self.top_diamond_alignment = self.alignments.len() - 1;
        }
    }
    pub fn find_deeparg_hit(&mut self, best_hit: String) {
        self.is_deeparg_hit = true;
        for alignment_index in 0..self.alignments.len(){
            if *self.alignments[alignment_index].get_name() == best_hit {
                self.top_deeparg_hit = Some(alignment_index);
                break;
            }
        }
    }
    pub fn update_amr_deeparg_hashmaps(&self, gene_map: &mut RcStringHashMap<usize>, 
            cdd_map: &mut RcStringHashMap<usize>, super_map: &mut RcStringHashMap<usize>,
            domain_name_refs: &mut StringHashMap<Rc<String>>) {
        if self.top_deeparg_hit.is_some() {
            self.alignments[self.top_deeparg_hit.unwrap()]
                .update_amr_hashmaps(gene_map, cdd_map, super_map, domain_name_refs);
        }
    }
    pub fn update_amr_restricted_diamond_hashmaps(&self, gene_map: &mut RcStringHashMap<usize>, 
            cdd_map: &mut RcStringHashMap<usize>, super_map: &mut RcStringHashMap<usize>,
            domain_name_refs: &mut StringHashMap<Rc<String>>) {
        if self.top_deeparg_hit.is_some() {
            self.alignments[self.top_diamond_alignment]
                .update_amr_hashmaps(gene_map, cdd_map, super_map, domain_name_refs);
        }
    }
    pub fn update_amr_diamond_hashmaps(&self, gene_map: &mut RcStringHashMap<usize>, 
            cdd_map: &mut RcStringHashMap<usize>, super_map: &mut RcStringHashMap<usize>,
            domain_name_refs: &mut StringHashMap<Rc<String>>) {
        self.alignments[self.top_diamond_alignment]
            .update_amr_hashmaps(gene_map, cdd_map, super_map, domain_name_refs);
    }
}
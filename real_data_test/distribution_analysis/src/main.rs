use std::collections::HashMap;
use std::hash::RandomState;
use std::cell::RefCell; 
use std::iter::zip;
use std::fs::File;
use std::io::Read;
use std::rc::Rc;

type RefVec = Vec<Rc<Reference>>;
type RefRc = Rc<Reference>;
type QueVec = Vec<RefCell<Query>>;
type QueRef = RefCell<Query>;
type StringHashMap<T> = HashMap<String, T, RandomState>;

static CDD_DIR: &str = "../../CDD_features/";
static DEEPARG_HIT_FILE: &str = "X.mapping.ARG";
static ALIGNMENT_FILE: &str = "X.align.daa.tsv";
static SAMPLE_ID_FILE: &str = "../real_samples.txt";
static PATH_TO_SS_OUTPUT: &str = "deeparg_results";
static PATH_TO_LS_OUTPUT: &str = "spades/deeparg_results";

fn main() {
    // Read content of real_samples.txt to find biosample ID
    let mut sample_id_file = File::open(SAMPLE_ID_FILE).unwrap();
    let mut sample_id_string = String::new();
    let _ = sample_id_file.read_to_string(&mut sample_id_string);
    drop(sample_id_file);
    let sample_id_vec: Vec<&str> = sample_id_string.split('\n').collect();

    // Look through CDD features and initialize a new Reference struct for each reference
    let mut ref_name_vec: Vec<String> = vec![];
    let mut gene_amr_hashmap: StringHashMap<usize> = HashMap::new();
    let mut cdd_amr_hashmap: StringHashMap<usize> = HashMap::new();
    let mut superfamily_amr_hashmap: StringHashMap<usize> = HashMap::new();
    let mut ref_struct_vec: RefVec = vec![];
    for part in 1..26usize {
        let mut cdd_info_file = File::open(format!("{CDD_DIR}Part{part}_hitdata.txt")).unwrap();
        let mut cdd_info_string = String::new();
        let _ = cdd_info_file.read_to_string(&mut cdd_info_string);
        let cdd_info_vec: Vec<&str> = cdd_info_string.split('\n').collect();
        for row_index in 1..cdd_info_vec.len() {
            let row: Vec<&str> = cdd_info_vec[row_index].split('\t').collect();
            let name: &str = row[0].split('>').last().unwrap();
            // cdd_acc.push(row[7].to_string());
            // superfamily_acc.push({
            //     if row[1] == "specific" {
            //         row[10].to_string()
            //     } else {
            //         row[7].to_string() 
            //     }}.to_string());
            if name == ref_name_vec.last().unwrap() {
                let mut ref_struct = Rc::try_unwrap(ref_struct_vec.pop().unwrap()).unwrap();
                ref_struct.add_to_domain(row);
                ref_struct_vec.push(Rc::new(ref_struct));
            }
            else {
                ref_name_vec.push(name.to_string());
                ref_struct_vec.push(Rc::new(Reference::new(row)));
            }
        }
    }
    
    // Make hashmap between reference name and reference struct to make matching constant time
    let ref_hashmap: StringHashMap<RefRc> = HashMap::from_iter(zip(ref_name_vec, ref_struct_vec));

    for sample_id in sample_id_vec {
        for identity in vec![30, 50, 80] {
            // Look through alignment output for short sequences
            let mut ss_query_name_vec: Vec<String> = vec![];
            let mut ss_query_struct_vec: QueVec = vec![];
            let mut full_ss_alignment_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_SS_OUTPUT, identity, ALIGNMENT_FILE)).unwrap();
            let mut full_ss_alignment_string = String::new();
            let _ = full_ss_alignment_file.read_to_string(&mut full_ss_alignment_string);
            drop(full_ss_alignment_file);
            let full_ss_alignment_vec: Vec<&str> = full_ss_alignment_string.split('\n').collect();
            for row_index in 0..full_ss_alignment_vec.len() {
                let row: Vec<&str> = full_ss_alignment_vec[row_index].split('\t').collect();
                if row[0] == ss_query_name_vec.last().unwrap() {
                    let mut query_struct = ss_query_struct_vec.last().unwrap().borrow_mut();
                    query_struct.add_alignment(row, &ref_hashmap);
                }
                else {
                    ss_query_name_vec.push(row[0].to_string());
                    ss_query_struct_vec.push(RefCell::new(Query::new(row, &ref_hashmap)));
                }
            }
            // Make hashmap between query name and query struct to make matching constant time
            let ss_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(
                zip(ss_query_name_vec, ss_query_struct_vec));

            // Look through deeparg hit file
            let mut full_ss_deeparg_hit_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_SS_OUTPUT, identity, DEEPARG_HIT_FILE)).unwrap();
            let mut full_ss_deeparg_hit_string = String::new();
            let _ = full_ss_deeparg_hit_file.read_to_string(&mut full_ss_deeparg_hit_string);
            drop(full_ss_deeparg_hit_file);
            let full_ss_deeparg_hit_vec: Vec<&str> = full_ss_deeparg_hit_string.split('\n').collect();
            for row_index in 1..full_ss_deeparg_hit_vec.len() {
                let row: Vec<&str> = full_ss_deeparg_hit_vec[row_index].split('\t').collect();
                let read_id = row[3].to_string();
                let best_hit = row[5].to_string();
                let mut query_struct = ss_query_hashmap.get(&read_id).unwrap().borrow_mut();
                query_struct.find_deeparg_hit(best_hit);
            }

            // Look through alignment output for long sequences
            let mut ls_query_name_vec: Vec<String> = vec![];
            let mut ls_query_struct_vec: QueVec = vec![];
            let mut full_ls_alignment_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_LS_OUTPUT, identity, ALIGNMENT_FILE)).unwrap();
            let mut full_ls_alignment_string = String::new();
            let _ = full_ls_alignment_file.read_to_string(&mut full_ls_alignment_string);
            drop(full_ls_alignment_file);
            let full_ls_alignment_vec: Vec<&str> = full_ls_alignment_string.split('\n').collect();
            for row_index in 0..full_ls_alignment_vec.len() {
                let row: Vec<&str> = full_ls_alignment_vec[row_index].split('\t').collect();
                if row[0] == ls_query_name_vec.last().unwrap() {
                    let mut query_struct = ls_query_struct_vec.last().unwrap().borrow_mut();
                    query_struct.add_alignment(row, &ref_hashmap);
                }
                else {
                    ls_query_name_vec.push(row[0].to_string());
                    ls_query_struct_vec.push(RefCell::new(Query::new(row, &ref_hashmap)));
                }
            }
            // Make hashmap between query name and query struct to make matching constant time
            let ls_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(zip(ls_query_name_vec, ls_query_struct_vec));

            // Look through deeparg hit file
            let mut full_ls_deeparg_hit_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_LS_OUTPUT, identity, DEEPARG_HIT_FILE)).unwrap();
            let mut full_ls_deeparg_hit_string = String::new();
            let _ = full_ls_deeparg_hit_file.read_to_string(&mut full_ls_deeparg_hit_string);
            drop(full_ls_deeparg_hit_file);
            let full_ls_deeparg_hit_vec: Vec<&str> = full_ls_deeparg_hit_string.split('\n').collect();
            for row_index in 1..full_ls_deeparg_hit_vec.len() {
                let row: Vec<&str> = full_ls_deeparg_hit_vec[row_index].split('\t').collect();
                let read_id = row[3].to_string();
                let best_hit = row[5].to_string();
                let mut query_struct = ls_query_hashmap.get(&read_id).unwrap().borrow_mut();
                query_struct.find_deeparg_hit(best_hit);
            }
        }
    }
}

#[derive(Debug)]
struct Domain {
    start: usize,
    end: usize,
    cdd_acc: String,
    superfamily_acc: String
}

#[derive(Debug)]
struct Reference {
    name: String,
    domains: Vec<Rc<Domain>>
}

impl Reference {
    fn new(row: Vec<&str>) -> Self {
        let name: &str = row[0].split('>').last().unwrap();
        Self {
            name: name.to_string(),
            domains: vec![Rc::new(Domain { 
                start: row[3].parse().unwrap(), 
                end: row[4].parse().unwrap(), 
                cdd_acc: row[7].to_string(), 
                superfamily_acc: if row[1] == "specific" {
                    row[10].to_string()
                } else {
                    row[7].to_string() 
                }
            })]
        }
    }
    fn add_to_domain(&mut self, row: Vec<&str>) {
        self.domains.push(Rc::new(Domain { 
            start: row[3].parse().unwrap(), 
            end: row[4].parse().unwrap(), 
            cdd_acc: row[7].to_string(), 
            superfamily_acc: if row[1] == "specific" {
                row[10].to_string()
            } else {
                row[7].to_string() 
            }
        }))
    }
    fn get_gene_name(&self) -> String {
        self.name.split("|").last().unwrap().to_string()
    }

    fn get_classification(&self) -> String {
        let fields: Vec<&str> = self.name.split("|").collect();
        fields[fields.len() - 2].to_string()
    }

    fn update_gene_amr_hashmap(&self, map: &mut StringHashMap<usize>){
        for domain in &self.domains {
            let key = format!("{}|{}|{}", 
                domain.cdd_acc, self.get_gene_name(), self.get_classification());
            if map.contains_key(&key) {
                *map.get_mut(&key).unwrap() += 1;
            }
        }
    }
}

#[derive(Debug)]
struct Alignment{
    matching_reference: Rc<Reference>,
    start: usize,
    end: usize,
    matching_domain: Vec<Rc<Domain>>
}

impl Alignment {
    fn new(row: Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) -> Self {
        let matching_reference = ref_hashmap.get(row[1]).unwrap().clone();
        let start = row[8].parse().unwrap();
        let end = row[9].parse().unwrap();
        let mut matching_domain = vec![];
        for domain in matching_reference.domains.iter() {
            if (domain.start >= start && domain.start < end) || (domain.end <= end && domain.end > start) {
                matching_domain.push(domain.clone())
            }
        }
        Self { matching_reference, start, end, matching_domain}
    }
    fn get_ref_name(&self) -> &String{
        &self.matching_reference.name
    }
}

#[derive(Debug)]
struct Query{
    name: String,
    alignments: Vec<Alignment>,
    top_diamond_alignment: usize,      //Should be the first alignment
    top_deeparg_hit: Option<usize>,    //Should be None until parsing DeepARG output
    is_deeparg_hit: bool
}

impl Query{
    fn new(row: Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) -> Self {
        Self { 
            name: row[0].to_string(), 
            alignments: vec![Alignment::new(row, ref_hashmap)],
            top_diamond_alignment: 0usize,
            top_deeparg_hit: None,
            is_deeparg_hit: false }
    }
    fn add_alignment(&mut self, row: Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) {
        self.alignments.push(Alignment::new(row, ref_hashmap));
    }
    fn find_deeparg_hit(&mut self, best_hit: String) {
        self.is_deeparg_hit = true;
        for alignment_index in 0..self.alignments.len(){
            if *self.alignments[alignment_index].get_ref_name() == best_hit {
                self.top_deeparg_hit = Some(alignment_index);
                break;
            }
        }
    }
}
use std::{fs::File, io::Read};
use boomphf::hashmap::BoomHashMap;

static CDD_DIR: &str = "../../CDD_features/";
static DEEPARG_MAP_FILE: &str = "X.mapping.ARG";
static ALIGNMENT_FILE: &str = "X.align.daa.tsv";
static SAMPLE_ID_FILE: &str = "../real_samples.txt";
static PATH_TO_DEEPARG_SS_OUTPUT: &str = "deeparg_results";
static PATH_TO_DEEPARG_LS_OUTPUT: &str = "spades/deeparg_results";

fn main() {
    // Read content of real_samples.txt to find biosample ID
    let mut sample_id_file = File::open(SAMPLE_ID_FILE).unwrap();
    let mut sample_id_string = String::new();
    let _ = sample_id_file.read_to_string(&mut sample_id_string);
    let sample_id_vec: Vec<&str> = sample_id_string.split('\n').collect();

    // Look through CDD features and initialize a new Reference struct for each reference
    let mut ref_name_vec: Vec<String> = vec![];
    let mut ref_struct_vec: Vec<Reference> = vec![];
    for part in 1..26usize {
        let mut cdd_info_file = File::open(format!("{CDD_DIR}Part{part}_hitdata.txt")).unwrap();
        let mut cdd_info_string = String::new();
        let _ = cdd_info_file.read_to_string(&mut cdd_info_string);
        let cdd_info_vec: Vec<&str> = cdd_info_string.split('\n').collect();
        for row_index in 1..cdd_info_vec.len() {
            let row: Vec<&str> = cdd_info_vec[row_index].split('\t').collect();
            let name: &str = row[0].split('>').last().unwrap();
            if name == ref_name_vec.last().unwrap() {
                let mut ref_struct = ref_struct_vec.pop().unwrap();
                ref_struct.add_to_domain(row);
                ref_struct_vec.push(ref_struct);
            }
            else {
                ref_name_vec.push(name.to_string());
                ref_struct_vec.push(Reference::new(row));
            }
        }
    }
    
    // Make hashmap between reference name and reference struct to make matching constant time
    let ref_hashmap: BoomHashMap<String, Reference> = BoomHashMap::new(ref_name_vec, ref_struct_vec);

    // Look through alignment output
    let mut query_name_vec: Vec<String> = vec![];


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
    domains: Vec<Domain>
}

impl Reference {
    fn new(row: Vec<&str>) -> Self {
        let name: &str = row[0].split('>').last().unwrap();
        Self {
            name: name.to_string(),
            domains: vec![Domain { 
                start: row[3].parse().unwrap(), 
                end: row[4].parse().unwrap(), 
                cdd_acc: row[7].to_string(), 
                superfamily_acc: if row[1] == "specific" {
                    row[10].to_string()
                } else {
                    row[7].to_string() 
                }
            }]
        }
    }
    fn add_to_domain(&mut self, row: Vec<&str>) {
        self.domains.push(Domain { 
            start: row[3].parse().unwrap(), 
            end: row[4].parse().unwrap(), 
            cdd_acc: row[7].to_string(), 
            superfamily_acc: if row[1] == "specific" {
                row[10].to_string()
            } else {
                row[7].to_string() 
            }
        });
    }
    fn get_gene_name(&self) -> String {
        self.name.split("|").last().unwrap().to_string()
    }

    fn get_classification(&self) -> String {
        let fields: Vec<&str> = self.name.split("|").collect();
        fields[fields.len() - 2].to_string()
    }
}

struct Alignment<'a>{
    matching_reference: &'a Reference,
    start: usize,
    end: usize,
    matching_domain: Vec<&'a Domain>
}

impl <'a> Alignment<'a> {
    fn new(row: Vec<&str>, ref_hashmap: &'a BoomHashMap<String, Reference>) -> Self {
        let matching_reference = ref_hashmap.get(row[1]).unwrap();
        let start = row[8].parse().unwrap();
        let end = row[9].parse().unwrap();
        let mut matching_domain = vec![];
        for domain in matching_reference.domains.iter() {
            if (domain.start >= start && domain.start < end) || (domain.end <= end && domain.end > start) {
                matching_domain.push(domain)
            }
        }
        Self { matching_reference, start, end, matching_domain}
    }
    
}

struct Query<'a>{
    name: String,
    alignments: Vec<&'a Alignment<'a>>
}

impl <'a> Query<'a> {
    fn new(row: Vec<&str>, ref_hashmap: &'a BoomHashMap<String, Reference>) -> Self {
        Self { 
            name: row[0].to_string(), 
            alignments: vec![] }
    }
}
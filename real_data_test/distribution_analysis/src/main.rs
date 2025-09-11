use std::collections::HashMap;
use std::hash::RandomState;
use std::io::{Read, Write};
use std::cell::RefCell; 
use std::iter::zip;
use std::fs::File;
use std::rc::Rc;

type RefVec = Vec<Rc<Reference>>;
type RefRc = Rc<Reference>;
type QueVec = Vec<RefCell<Query>>;
type QueRef = RefCell<Query>;
type StringHashMap<T> = HashMap<String, T, RandomState>;
type RcStringHashMap<T> = HashMap<Rc<String>, T, RandomState>;

static CDD_DIR: &str = "../../CDD_features/";
static REF_LOC: &str = "../../data/database/v2/features.fasta";
static DEEPARG_HIT_FILE: &str = "X.mapping.ARG";
static ALIGNMENT_FILE: &str = "X.align.daa.tsv";
static SAMPLE_ID_FILE: &str = "../real_samples.txt";
static PATH_TO_SS: &str = "deeparg_results";
static PATH_TO_LS: &str = "spades/deeparg_results";

static COMPARISON_OUTPUT: &str = "comparison.txt";

fn main() {
    // Read content of real_samples.txt to find biosample ID
    let mut sample_id_file = File::open(SAMPLE_ID_FILE).unwrap();
    let mut sample_id_string = String::new();
    let _ = sample_id_file.read_to_string(&mut sample_id_string);
    drop(sample_id_file);
    let sample_id_vec: Vec<&str> = sample_id_string.split('\n').collect();

    // Make distribution hashmap for entire database
    let mut domain_name_refs: StringHashMap<Rc<String>> = HashMap::new();
    let mut db_gene_amr_hashmap: RcStringHashMap<f64> = HashMap::new();
    let mut db_cdd_amr_hashmap: RcStringHashMap<f64> = HashMap::new();
    let mut db_super_amr_hashmap: RcStringHashMap<f64> = HashMap::new();

    // Look through CDD features and initialize a new Reference struct for each reference
    let mut ref_name_vec: Vec<String> = vec![];
    let mut ref_struct_vec: RefVec = vec![];
    for part in 1..26usize {
        let mut cdd_info_file = File::open(format!("{CDD_DIR}Part{part}_hitdata.txt")).unwrap();
        let mut cdd_info_string = String::new();
        let _ = cdd_info_file.read_to_string(&mut cdd_info_string);
        let cdd_info_vec: Vec<&str> = cdd_info_string.split('\n').collect();
        for row_index in 1..cdd_info_vec.len()-1 {
            let row: Vec<&str> = cdd_info_vec[row_index].split('\t').collect();
            let name: &str = row[0].split('>').last().unwrap();
            if ref_name_vec.len() > 0 && name == ref_name_vec.last().unwrap() {
                let mut ref_struct = Rc::try_unwrap(ref_struct_vec.pop().unwrap()).unwrap();
                ref_struct.add_to_domain(row);
                ref_struct_vec.push(Rc::new(ref_struct));
            }
            else {
                if part > 1 || row_index > 1 {
                    ref_struct_vec.last().unwrap().update_amr_hashmaps(
                        &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, 
                        &mut db_super_amr_hashmap, &mut domain_name_refs);
                }
                ref_name_vec.push(name.to_string());
                ref_struct_vec.push(Rc::new(Reference::new(row)));
            }
        }
    }

    ref_struct_vec.last().unwrap().update_amr_hashmaps(
        &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, 
        &mut db_super_amr_hashmap, &mut domain_name_refs);
    let mut ref_file = File::open(REF_LOC).unwrap();
    let mut ref_string = String::new();
    let _ = ref_file.read_to_string(&mut ref_string);
    drop(ref_file);
    let ref_vec: Vec<&str> = ref_string.split('\n').collect();
    for ref_index in 0..(ref_vec.len()-1)/2 {
        let ref_name = &ref_vec[ref_index*2][1..];
        let ref_seq = ref_vec[ref_index*2+1];
        if !ref_name_vec.contains(&ref_name.to_string()) {
            ref_name_vec.push(ref_name.to_string());
            ref_struct_vec.push(
                Rc::new(Reference {
                    name: ref_name.to_string(),
                    domains: vec![Rc::new(Domain { 
                        start: 0usize, 
                        end: ref_seq.len()-1, 
                        cdd_acc: ref_name.split('|').last().unwrap().to_string(), 
                        super_acc: ref_name.split('|').last().unwrap().to_string()
                    })]
                })
            );
            ref_struct_vec.last().unwrap().update_amr_hashmaps(
                &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, 
                &mut db_super_amr_hashmap, &mut domain_name_refs);
        }
    }

    // Instances when deeparg hit != diamond alignment. Tuple has sample; identity; 
    // model type; diamond's gene_amr info, cdd_amr info, and super_amr info; deeparg's 
    // gene_amr info, cdd_amr info, and super_amr info where info is key, distr
    // in ref, distr in diamond, distr in diamond(restricted) and distr in deeparg
    let mut deeparg_not_equal_diamond: HashMap<(&str, i32, &str, 
        [String; 6]), usize, RandomState> = HashMap::new();

    // Make hashmap between reference name and reference struct to make matching constant time
    let ref_hashmap: StringHashMap<RefRc> = HashMap::from_iter(zip(ref_name_vec, ref_struct_vec));

    for sample_id in sample_id_vec {
    println!("{}", sample_id);
        for identity in vec![30, 50, 80] {
            println!("{}", identity);
            println!("SS");
            // Make distribution hashmap for analysis 3
            let mut gene_amr_hashmap3: RcStringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap3: RcStringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap3: RcStringHashMap<f64> = HashMap::new();

            // Look through alignment output for short sequences
            let mut ss_query_name_vec: Vec<String> = vec![];
            let mut ss_query_struct_vec: QueVec = vec![];
            let mut full_ss_alignment_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_SS, identity, ALIGNMENT_FILE)).unwrap();
            let mut full_ss_alignment_string = String::new();
            let _ = full_ss_alignment_file.read_to_string(&mut full_ss_alignment_string);
            drop(full_ss_alignment_file);
            let full_ss_alignment_vec: Vec<&str> = full_ss_alignment_string.split('\n').collect();
            for row_index in 0..full_ss_alignment_vec.len()-1 {
                let row: Vec<&str> = full_ss_alignment_vec[row_index].split('\t').collect();
                if ss_query_name_vec.len() > 0 && row[0] == ss_query_name_vec.last().unwrap() {
                    let mut query_struct = ss_query_struct_vec.last().unwrap().borrow_mut();
                    query_struct.add_alignment(&row, &ref_hashmap);
                }
                else {
                    ss_query_name_vec.push(row[0].to_string());
                    ss_query_struct_vec.push(RefCell::new(Query::new(&row, &ref_hashmap)));
                    ss_query_struct_vec.last().unwrap().borrow().update_amr_hashmaps_v3(
                        &mut gene_amr_hashmap3, &mut cdd_amr_hashmap3, 
                        &mut super_amr_hashmap3, &mut domain_name_refs);
                }
            }
            // Make hashmap between query name and query struct to make matching constant time
            let ss_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(
                zip(ss_query_name_vec, ss_query_struct_vec));

            // Make distribution hashmap for analysis 1
            let mut gene_amr_hashmap1: RcStringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap1: RcStringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap1: RcStringHashMap<f64> = HashMap::new();

            // Make distribution hashmap for analysis 2
            let mut gene_amr_hashmap2: RcStringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap2: RcStringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap2: RcStringHashMap<f64> = HashMap::new();

            // Look through deeparg hit file
            let mut full_ss_deeparg_hit_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_SS, identity, DEEPARG_HIT_FILE)).unwrap();
            let mut full_ss_deeparg_hit_string = String::new();
            let _ = full_ss_deeparg_hit_file.read_to_string(&mut full_ss_deeparg_hit_string);
            drop(full_ss_deeparg_hit_file);
            let full_ss_deeparg_hit_vec: Vec<&str> = full_ss_deeparg_hit_string.split('\n').collect();
            for row_index in 1..full_ss_deeparg_hit_vec.len()-1 {
                let row: Vec<&str> = full_ss_deeparg_hit_vec[row_index].split('\t').collect();
                let read_id = row[3].to_string();
                let best_hit = row[5].to_string();
                let mut query_struct = ss_query_hashmap.get(&read_id).unwrap().borrow_mut();
                query_struct.find_deeparg_hit(best_hit);
                query_struct.update_amr_hashmaps_v1(
                    &mut gene_amr_hashmap1, &mut cdd_amr_hashmap1, 
                    &mut super_amr_hashmap1, &mut domain_name_refs);
                query_struct.update_amr_hashmaps_v2(
                    &mut gene_amr_hashmap2, &mut cdd_amr_hashmap2, 
                    &mut super_amr_hashmap2, &mut domain_name_refs);
            }

            for query in ss_query_hashmap.values() {
                let borrowed = query.borrow();
                if !borrowed.is_deeparg_hit {
                    continue;
                }
                if borrowed.top_diamond_alignment == borrowed.top_deeparg_hit.unwrap() {
                    continue;
                }

                let dia_doms = borrowed.get_top_diamond_alignment_domain_identifiers();
                let dee_doms = borrowed.get_top_deeparg_hit_domain_identifiers().unwrap();
                for (dia_gene_amr, dia_cdd_amr, dia_super_amr) in &dia_doms {
                    for (dee_gene_amr, dee_cdd_amr, dee_super_amr) in &dee_doms {
                        if *dia_gene_amr == *dee_gene_amr {
                            continue;
                        }
                        if dia_doms.contains(&(dee_gene_amr.clone(), dee_cdd_amr.clone(), dee_super_amr.clone())) &&
                            dee_doms.contains(&(dia_gene_amr.clone(), dia_cdd_amr.clone(), dia_super_amr.clone())) {
                            continue;
                        }
                        let dee_gene_amr_info = 
                            get_domain_info(dee_gene_amr.to_string(), &gene_amr_hashmap1,
                                &gene_amr_hashmap2, &gene_amr_hashmap3, 
                                &db_gene_amr_hashmap, &mut domain_name_refs);
                        let dee_cdd_amr_info = 
                            get_domain_info(dee_cdd_amr.to_string(), &cdd_amr_hashmap1,
                                &cdd_amr_hashmap2, &cdd_amr_hashmap3, 
                                &db_cdd_amr_hashmap, &mut domain_name_refs);
                        let dee_super_amr_info = 
                            get_domain_info(dee_super_amr.to_string(), &super_amr_hashmap1,
                                &super_amr_hashmap2, &super_amr_hashmap3, 
                                &db_super_amr_hashmap, &mut domain_name_refs);
                        let dia_gene_amr_info = 
                            get_domain_info(dia_gene_amr.to_string(), &gene_amr_hashmap1,
                                &gene_amr_hashmap2, &gene_amr_hashmap3, 
                                &db_gene_amr_hashmap, &mut domain_name_refs);
                        let dia_cdd_amr_info = 
                            get_domain_info(dia_cdd_amr.to_string(), &cdd_amr_hashmap1,
                                &cdd_amr_hashmap2, &cdd_amr_hashmap3, 
                                &db_cdd_amr_hashmap, &mut domain_name_refs);
                        let dia_super_amr_info = 
                            get_domain_info(dia_super_amr.to_string(), &super_amr_hashmap1,
                                &super_amr_hashmap2, &super_amr_hashmap3, 
                                &db_super_amr_hashmap, &mut domain_name_refs);
                        let key = (sample_id, identity, "SS", [
                            dia_gene_amr_info, dia_cdd_amr_info, dia_super_amr_info,
                            dee_gene_amr_info, dee_cdd_amr_info, dee_super_amr_info,
                        ]);
                        if !deeparg_not_equal_diamond.contains_key(&key) {
                            deeparg_not_equal_diamond.insert(key, 1);
                        }
                        else {
                            *deeparg_not_equal_diamond.get_mut(&key).unwrap() += 1;
                        }
                    }
                }
            }

            println!("LS");
            // Make distribution hashmap for analysis 3
            let mut gene_amr_hashmap3: RcStringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap3: RcStringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap3: RcStringHashMap<f64> = HashMap::new();

            // Look through alignment output for long sequences
            let mut ls_query_name_vec: Vec<String> = vec![];
            let mut ls_query_struct_vec: QueVec = vec![];
            let mut full_ls_alignment_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_LS, identity, ALIGNMENT_FILE)).unwrap();
            let mut full_ls_alignment_string = String::new();
            let _ = full_ls_alignment_file.read_to_string(&mut full_ls_alignment_string);
            drop(full_ls_alignment_file);
            let full_ls_alignment_vec: Vec<&str> = full_ls_alignment_string.split('\n').collect();
            for row_index in 0..full_ls_alignment_vec.len()-1 {
                let row: Vec<&str> = full_ls_alignment_vec[row_index].split('\t').collect();
                if ls_query_name_vec.len() > 0 && row[0] == ls_query_name_vec.last().unwrap() {
                    let mut query_struct = ls_query_struct_vec.last().unwrap().borrow_mut();
                    query_struct.add_alignment(&row, &ref_hashmap);
                }
                else {
                    ls_query_name_vec.push(row[0].to_string());
                    ls_query_struct_vec.push(RefCell::new(Query::new(&row, &ref_hashmap)));
                    ls_query_struct_vec.last().unwrap().borrow().update_amr_hashmaps_v3(
                        &mut gene_amr_hashmap3, &mut cdd_amr_hashmap3, 
                        &mut super_amr_hashmap3, &mut domain_name_refs);
                }
            }
            // Make hashmap between query name and query struct to make matching constant time
            let ls_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(zip(ls_query_name_vec, ls_query_struct_vec));

            // Make distribution hashmap for analysis 1
            let mut gene_amr_hashmap1: RcStringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap1: RcStringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap1: RcStringHashMap<f64> = HashMap::new();

            // Make distribution hashmap for analysis 2
            let mut gene_amr_hashmap2: RcStringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap2: RcStringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap2: RcStringHashMap<f64> = HashMap::new();

            // Look through deeparg hit file
            let mut full_ls_deeparg_hit_file = File::open(format!(
                "../samples/{}/{}/arg_alignment_identity_{}/{}",
                sample_id, PATH_TO_LS, identity, DEEPARG_HIT_FILE)).unwrap();
            let mut full_ls_deeparg_hit_string = String::new();
            let _ = full_ls_deeparg_hit_file.read_to_string(&mut full_ls_deeparg_hit_string);
            drop(full_ls_deeparg_hit_file);
            let full_ls_deeparg_hit_vec: Vec<&str> = full_ls_deeparg_hit_string.split('\n').collect();
            for row_index in 1..full_ls_deeparg_hit_vec.len()-1 {
                let row: Vec<&str> = full_ls_deeparg_hit_vec[row_index].split('\t').collect();
                let read_id = row[3].to_string();
                let best_hit = row[5].to_string();

                // Mark query struct that matches with read_id as deeparg hit
                // And save alignment that matches with best hit in main_deeparg_hit
                let mut query_struct = ls_query_hashmap.get(&read_id).unwrap().borrow_mut();
                query_struct.find_deeparg_hit(best_hit);
                query_struct.update_amr_hashmaps_v1(
                    &mut gene_amr_hashmap1, &mut cdd_amr_hashmap1, 
                    &mut super_amr_hashmap1, &mut domain_name_refs);
                query_struct.update_amr_hashmaps_v2(
                    &mut gene_amr_hashmap2, &mut cdd_amr_hashmap2, 
                    &mut super_amr_hashmap2, &mut domain_name_refs);
            }

            for query in ls_query_hashmap.values() {
                let borrowed = query.borrow();
                if !borrowed.is_deeparg_hit {
                    continue;
                }
                if borrowed.top_diamond_alignment == borrowed.top_deeparg_hit.unwrap() {
                    continue;
                }

                let dia_doms = borrowed.get_top_diamond_alignment_domain_identifiers();
                let dee_doms = borrowed.get_top_deeparg_hit_domain_identifiers().unwrap();
                for (dia_gene_amr, dia_cdd_amr, dia_super_amr) in &dia_doms {
                    for (dee_gene_amr, dee_cdd_amr, dee_super_amr) in &dee_doms {
                        if *dia_gene_amr == *dee_gene_amr {
                            continue;
                        }
                        if dia_doms.contains(&(dee_gene_amr.clone(), dee_cdd_amr.clone(), dee_super_amr.clone())) &&
                            dee_doms.contains(&(dia_gene_amr.clone(), dia_cdd_amr.clone(), dia_super_amr.clone())) {
                            continue;
                        }
                        let dee_gene_amr_info = 
                            get_domain_info(dee_gene_amr.to_string(), &gene_amr_hashmap1,
                                &gene_amr_hashmap2, &gene_amr_hashmap3, 
                                &db_gene_amr_hashmap, &mut domain_name_refs);
                        let dee_cdd_amr_info = 
                            get_domain_info(dee_cdd_amr.to_string(), &cdd_amr_hashmap1,
                                &cdd_amr_hashmap2, &cdd_amr_hashmap3, 
                                &db_cdd_amr_hashmap, &mut domain_name_refs);
                        let dee_super_amr_info = 
                            get_domain_info(dee_super_amr.to_string(), &super_amr_hashmap1,
                                &super_amr_hashmap2, &super_amr_hashmap3, 
                                &db_super_amr_hashmap, &mut domain_name_refs);
                        let dia_gene_amr_info = 
                            get_domain_info(dia_gene_amr.to_string(), &gene_amr_hashmap1,
                                &gene_amr_hashmap2, &gene_amr_hashmap3, 
                                &db_gene_amr_hashmap, &mut domain_name_refs);
                        let dia_cdd_amr_info = 
                            get_domain_info(dia_cdd_amr.to_string(), &cdd_amr_hashmap1,
                                &cdd_amr_hashmap2, &cdd_amr_hashmap3, 
                                &db_cdd_amr_hashmap, &mut domain_name_refs);
                        let dia_super_amr_info = 
                            get_domain_info(dia_super_amr.to_string(), &super_amr_hashmap1,
                                &super_amr_hashmap2, &super_amr_hashmap3, 
                                &db_super_amr_hashmap, &mut domain_name_refs);
                        let key = (sample_id, identity, "LS", [
                            dia_gene_amr_info, dia_cdd_amr_info, dia_super_amr_info,
                            dee_gene_amr_info, dee_cdd_amr_info, dee_super_amr_info,
                        ]);
                        if !deeparg_not_equal_diamond.contains_key(&key) {
                            deeparg_not_equal_diamond.insert(key, 1);
                        }
                        else {
                            *deeparg_not_equal_diamond.get_mut(&key).unwrap() += 1;
                        }
                    }
                }
            }
        }
    }

    // Create output files
    let mut output = File::create(COMPARISON_OUTPUT).unwrap();
    //Tuple has sample; identity; 
    // model type; diamond's gene_amr info, cdd_amr info, and super_amr info; deeparg's 
    // gene_amr info, cdd_amr info, and super_amr info where info is key, distr
    // in ref, distr in diamond, distr in diamond(restricted) and distr in deeparg
    let _ = output.write(b"Sample\tAlignment Identity\tModel\tDiamond alignment cdd|gene|amr key\t");
    let _ = output.write(b"Diamond alignment cdd|gene|amr reference distribution\t");
    let _ = output.write(b"Diamond alignment cdd|gene|amr diamond distribution\t");
    let _ = output.write(b"Diamond alignment cdd|gene|amr diamond restricted distribution\t");
    let _ = output.write(b"Diamond alignment cdd|gene|amr deeparg distribution\t");
    let _ = output.write(b"Diamond alignment cdd|amr key\t");
    let _ = output.write(b"Diamond alignment cdd|amr reference distribution\t");
    let _ = output.write(b"Diamond alignment cdd|amr diamond distribution\t");
    let _ = output.write(b"Diamond alignment cdd|amr diamond restricted distribution\t");
    let _ = output.write(b"Diamond alignment cdd|amr deeparg distribution\t");
    let _ = output.write(b"Diamond alignment super|amr key\t");
    let _ = output.write(b"Diamond alignment super|amr reference distribution\t");
    let _ = output.write(b"Diamond alignment super|amr diamond distribution\t");
    let _ = output.write(b"Diamond alignment super|amr diamond restricted distribution\t");
    let _ = output.write(b"Diamond alignment super|amr deeparg distribution\t");
    let _ = output.write(b"DeepARG hit cdd|gene|amr key\t");
    let _ = output.write(b"DeepARG hit cdd|gene|amr reference distribution\t");
    let _ = output.write(b"DeepARG hit cdd|gene|amr diamond distribution\t");
    let _ = output.write(b"DeepARG hit cdd|gene|amr diamond restricted distribution\t");
    let _ = output.write(b"DeepARG hit cdd|gene|amr deeparg distribution\t");
    let _ = output.write(b"DeepARG hit cdd|amr key\t");
    let _ = output.write(b"DeepARG hit cdd|amr reference distribution\t");
    let _ = output.write(b"DeepARG hit cdd|amr diamond distribution\t");
    let _ = output.write(b"DeepARG hit cdd|amr diamond restricted distribution\t");
    let _ = output.write(b"DeepARG hit cdd|amr deeparg distribution\t");
    let _ = output.write(b"DeepARG hit super|amr key\t");
    let _ = output.write(b"DeepARG hit super|amr reference distribution\t");
    let _ = output.write(b"DeepARG hit super|amr diamond distribution\t");
    let _ = output.write(b"DeepARG hit super|amr diamond restricted distribution\t");
    let _ = output.write(b"DeepARG hit super|amr deeparg distribution\t");
    let _ = output.write(b"Diamond alignment cdd|gene|amr/DeepARG hit cdd|gene|amr deeparg count\n");
    write_hashmap_to_file(&mut output, deeparg_not_equal_diamond);
}

fn get_domain_info(dom: String, 
        hashmap_1: &RcStringHashMap<f64>, hashmap_2: &RcStringHashMap<f64>, 
        hashmap_3: &RcStringHashMap<f64>, hashmap_db: &RcStringHashMap<f64>, 
        domain_name_refs: &mut StringHashMap<Rc<String>>) -> String {
    let dom_refs = domain_name_refs.get(&dom).unwrap();
    if hashmap_2.get(dom_refs).is_none() && hashmap_1.get(dom_refs).is_none() {
        panic!("Messed up with {dom}")
    }
    format!("{}\t{}\t{}\t{}\t{}", dom.clone(), *hashmap_db.get(dom_refs).unwrap(), 
        *hashmap_3.get(dom_refs).unwrap_or(&0.0), *hashmap_2.get(dom_refs).unwrap_or(&0.0),
        *hashmap_1.get(dom_refs).unwrap_or(&0.0))
}

fn write_hashmap_to_file(file: &mut File, map: HashMap<(&str, i32, &str, 
        [String; 6]), usize, RandomState>) {
    for row in map.iter() {
        let _ = file.write(format!("{}\t{}\t{}\t", (*row.0).0,
            (*row.0).1, (*row.0).2).as_bytes());
        for val in &(*row.0).3 {
            let _ = file.write(format!("{}\t", val).as_bytes());
        }
        let _ = file.write(format!("{}\n", row.1).as_bytes());
    }
}

#[derive(Debug)]
struct Domain {
    start: usize,
    end: usize,
    cdd_acc: String,
    super_acc: String
}

trait DomainContainer{
    // We want a length-weighted distribution 
    fn get_domains(&self) -> &Vec<Rc<Domain>>;
    fn get_name(&self) -> &String;
    fn get_restrictions(&self) -> Option<(usize, usize)>;

    fn get_calculated_weights(&self) -> Vec<Vec<f64>> {
        // Some domains may overlap, so we'll divide the weight 
        // of the overlapping region by the degree of overlap
        let mut domains_weight_dividend: Vec<Vec<f64>> = vec![];
        for index in 0..self.get_domains().len() {
            let curr_domain = self.get_domains()[index].clone();
            let curr_start = curr_domain.start;
            let curr_end = curr_domain.end;
            let curr_length = curr_end - curr_start + 1;
            drop(curr_domain);
            let mut curr_domain_weight_dividend = vec![1.0; curr_end-curr_start + 1];
            if index == 0 {
                domains_weight_dividend.push(curr_domain_weight_dividend);
                continue;
            }
            for prev_index in 0..index {
                let prev_domain = self.get_domains()[prev_index].clone();
                let prev_start = prev_domain.start;
                let prev_end = prev_domain.end;
                let prev_length = prev_end - prev_start + 1;
                drop(prev_domain);
                if prev_start >= curr_start && prev_end <= curr_end {
                    // prev domain is completely within curr domain
                    domains_weight_dividend[prev_index] = domains_weight_dividend[prev_index]
                        .iter().map(|x| 1.0 + x).collect();
                    let curr_overlap_start = prev_start - curr_start;
                    let curr_overlap_end = prev_end - curr_start + 1;
                    for pos in curr_overlap_start..curr_overlap_end {
                        curr_domain_weight_dividend[pos] += 1.0;
                    }
                }
                else if prev_start >= curr_start && prev_start <= curr_end{
                    // the start of prev overlaps with the end of curr
                    let curr_overlap_start = prev_start - curr_start;
                    let prev_overlap_end = curr_end - prev_start + 1;
                    for pos in curr_overlap_start..curr_length {
                        curr_domain_weight_dividend[pos] += 1.0;
                    }
                    for pos in 0..prev_overlap_end {
                        domains_weight_dividend[prev_index][pos] += 1.0;
                    }
                }

                else if curr_start >= prev_start && curr_end <= prev_end {
                    // curr domain is completely within prev domain
                    curr_domain_weight_dividend = curr_domain_weight_dividend
                        .iter().map(|x| 1.0 + x).collect();
                    let prev_overlap_start = curr_start - prev_start;
                    let prev_overlap_end = curr_end - prev_start + 1;
                    for pos in prev_overlap_start..prev_overlap_end {
                        domains_weight_dividend[prev_index][pos] += 1.0;
                    }
                }
                
                else if curr_start >= prev_start && curr_start <= prev_end {
                    // the start of curr overlaps with the end of prev
                    let prev_overlap_start = curr_start - prev_start;
                    let curr_overlap_end = prev_end - curr_start + 1;
                    for pos in prev_overlap_start..prev_length {
                        domains_weight_dividend[prev_index][pos] += 1.0;
                    }
                    for pos in 0..curr_overlap_end {
                        curr_domain_weight_dividend[pos] += 1.0;
                    }
                }
            }
            domains_weight_dividend.push(curr_domain_weight_dividend);
        }
        let mut domains_weight: Vec<Vec<f64>> = vec![];
        for index in 0..self.get_domains().len() {
            domains_weight.push(domains_weight_dividend[index].iter()
                .map(|x| 1.0/x).collect()
            );
        }
        domains_weight
    }

    fn get_gene_name(&self) -> String {
        self.get_name().split("|").last().unwrap().to_string()
    }

    fn get_classification(&self) -> String {
        let fields: Vec<&str> = self.get_name().split("|").collect();
        fields[fields.len() - 2].to_string()
    }

    fn get_domain_identifiers(&self) -> Vec<(String, String, String)> {
        let mut domain_identifiers: Vec<(String, String, String)> = vec![];
        for domain in self.get_domains() {
            let gene_key = format!("{}|{}|{}", 
                domain.cdd_acc, self.get_gene_name().to_uppercase(), self.get_classification());
            let cdd_key = format!("{}|{}", 
                domain.cdd_acc, self.get_classification());
            let super_key = format!("{}|{}", 
                domain.super_acc, self.get_classification());
            domain_identifiers.push((gene_key, cdd_key, super_key));
        }
        return domain_identifiers
    }

    fn update_amr_hashmaps(&self, gene_map: &mut RcStringHashMap<f64>, 
            cdd_map: &mut RcStringHashMap<f64>, super_map: &mut RcStringHashMap<f64>,
            domain_name_refs: &mut StringHashMap<Rc<String>>){
        let domains_weight = self.get_calculated_weights();
        for domain_index in 0..self.get_domains().len() {
            let domain = self.get_domains()[domain_index].clone();
            let weight = if self.get_restrictions().is_some() {
                let restrict_start = self.get_restrictions().unwrap().0;
                let restrict_end = self.get_restrictions().unwrap().1;
                if restrict_start >= domain.start && restrict_end <= domain.end {
                    domains_weight[domain_index][restrict_start-domain.start..
                        restrict_end-domain.start+1].iter().sum::<f64>()
                }
                else if restrict_start >= domain.start && restrict_start <= domain.end{
                    domains_weight[domain_index][restrict_start-domain.start..
                        ].iter().sum::<f64>()
                }
                else if restrict_end >= domain.start && restrict_end <= domain.end{
                    domains_weight[domain_index][..
                        restrict_end-domain.start+1].iter().sum::<f64>()
                }
                else if restrict_start < domain.start && restrict_end > domain.end{
                    domains_weight[domain_index].iter().sum::<f64>()
                }
                else {
                    0.0
                }
            }
            else {
                domains_weight[domain_index].iter().sum::<f64>()
            };
            let gene_key = format!("{}|{}|{}", 
                domain.cdd_acc, self.get_gene_name().to_uppercase(), self.get_classification());
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
                *gene_map.get_mut(&gene_key_ref).unwrap() += weight
                
            }
            else {
                gene_map.insert(gene_key_ref, weight);
            }
            let cdd_key = format!("{}|{}", 
                domain.cdd_acc, self.get_classification());
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
                *cdd_map.get_mut(&cdd_key_ref).unwrap() += weight;
            }
            else {
                cdd_map.insert(cdd_key_ref, weight);
            }
            let super_key = format!("{}|{}", 
                domain.super_acc, self.get_classification());
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
                *super_map.get_mut(&super_key_ref).unwrap() += weight;
            }
            else {
                super_map.insert(super_key_ref, weight);
            }
        }
    }
}

#[derive(Debug)]
struct Reference {
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
    fn get_restrictions(&self) -> Option<(usize, usize)> {
        None
    }
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
                super_acc: if row[1] == "specific" {
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
            super_acc: if row[1] == "specific" {
                row[10].to_string()
            } else {
                row[7].to_string() 
            }
        }))
    }
    
}

#[derive(Debug)]
struct Alignment{
    matching_reference: Rc<Reference>,
    start: usize,
    end: usize,
    bitscore: f64,
    domains: Vec<Rc<Domain>>
}

impl DomainContainer for Alignment {
    fn get_domains(&self) -> &Vec<Rc<Domain>> {
        &self.domains
    }
    fn get_name(&self) -> &String {
        &self.matching_reference.name
    }
    fn get_restrictions(&self) -> Option<(usize, usize)> {
        Some((self.start, self.end))
    }
}

impl Alignment {
    fn new(row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) -> Self {
        let matching_reference = ref_hashmap.get(row[1]).unwrap().clone();
        let start = row[8].parse().unwrap();
        let end = row[9].parse().unwrap();
        let bitscore = row[11].parse().unwrap();
        let mut domains = vec![];
        for domain in matching_reference.domains.iter() {
            if (domain.start >= start && domain.start < end) || (domain.end <= end && domain.end > start) {
                domains.push(domain.clone())
            }
        }
        Self { matching_reference, start, end, bitscore, domains}
    }
    fn get_ref_name(&self) -> &String{
        &self.matching_reference.name
    }
}

#[derive(Debug)]
struct Query{
    alignments: Vec<Alignment>,
    top_diamond_alignment: usize,      //Should be the first alignment
    top_deeparg_hit: Option<usize>,    //Should be None until parsing DeepARG output
    is_deeparg_hit: bool
}

impl Query{
    fn new(row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) -> Self {
        Self { 
            alignments: vec![Alignment::new(row, ref_hashmap)],
            top_diamond_alignment: 0usize,
            top_deeparg_hit: None,
            is_deeparg_hit: false }
    }

    fn get_top_diamond_alignment_domain_identifiers(&self) -> Vec<(String, String, String)> {
        self.alignments[self.top_diamond_alignment].get_domain_identifiers()
    }

    fn get_top_deeparg_hit_domain_identifiers(&self) -> Result<Vec<(String, String, String)>, &str> {
        if self.top_deeparg_hit.is_none() {
            Err("Query is not part of DeepARG output")
        }
        else {
            Ok(self.alignments[self.top_deeparg_hit.unwrap()].get_domain_identifiers())
        }
    }

    fn add_alignment(&mut self, row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) {
        self.alignments.push(Alignment::new(row, ref_hashmap));
        if self.alignments.last().unwrap().bitscore > self.alignments[self.top_diamond_alignment].bitscore {
            self.top_diamond_alignment = self.alignments.len() - 1;
        }
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
    fn update_amr_hashmaps_v1(&self, gene_map: &mut RcStringHashMap<f64>, 
            cdd_map: &mut RcStringHashMap<f64>, super_map: &mut RcStringHashMap<f64>,
            domain_name_refs: &mut StringHashMap<Rc<String>>) {
        if self.top_deeparg_hit.is_some() {
            self.alignments[self.top_deeparg_hit.unwrap()]
                .update_amr_hashmaps(gene_map, cdd_map, super_map, domain_name_refs);
        }
    }
    fn update_amr_hashmaps_v2(&self, gene_map: &mut RcStringHashMap<f64>, 
            cdd_map: &mut RcStringHashMap<f64>, super_map: &mut RcStringHashMap<f64>,
            domain_name_refs: &mut StringHashMap<Rc<String>>) {
        if self.top_deeparg_hit.is_some() {
            self.alignments[self.top_diamond_alignment]
                .update_amr_hashmaps(gene_map, cdd_map, super_map, domain_name_refs);
        }
    }
    fn update_amr_hashmaps_v3(&self, gene_map: &mut RcStringHashMap<f64>, 
            cdd_map: &mut RcStringHashMap<f64>, super_map: &mut RcStringHashMap<f64>,
            domain_name_refs: &mut StringHashMap<Rc<String>>) {
        self.alignments[self.top_diamond_alignment]
            .update_amr_hashmaps(gene_map, cdd_map, super_map, domain_name_refs);
    }
}
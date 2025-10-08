use std::collections::HashMap;
use std::hash::RandomState;
use std::io::{Read, Write};
use std::cell::RefCell; 
use std::iter::zip;
use std::fs::File;
use std::rc::Rc;

use domain::DomainContainer;
use reference::Reference;
use domain::Domain;
use query::Query;

pub mod reference;
pub mod alignment;
pub mod domain;
pub mod query;

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
static AMR_SWITCH_DISTRIBUTION: &str = "amr_switch_distribution.txt";

fn main() {
    // Read content of real_samples.txt to find biosample ID
    let mut sample_id_file = File::open(SAMPLE_ID_FILE).unwrap();
    let mut sample_id_string = String::new();
    let _ = sample_id_file.read_to_string(&mut sample_id_string);
    drop(sample_id_file);
    let sample_id_vec: Vec<&str> = sample_id_string.split('\n').collect();

    // Make distribution hashmap for entire database
    let mut domain_name_refs: StringHashMap<Rc<String>> = HashMap::new();
    let mut db_gene_amr_hashmap: RcStringHashMap<usize> = HashMap::new();
    let mut db_cdd_amr_hashmap: RcStringHashMap<usize> = HashMap::new();
    let mut db_super_amr_hashmap: RcStringHashMap<usize> = HashMap::new();

    // Count amr classes in database
    let mut db_amr_count: StringHashMap<usize> = HashMap::new();

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
                let amr = ref_struct_vec.last().unwrap().get_classification();
                if db_amr_count.contains_key(&amr) {
                    *db_amr_count.get_mut(&amr).unwrap() += 1;
                }
                else {
                    db_amr_count.insert(amr, 1);
                }
            }
        }
    }

    ref_struct_vec.last().unwrap().update_amr_hashmaps(
        &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, 
        &mut db_super_amr_hashmap, &mut domain_name_refs);

    // Add additional references that don't have domains
    let mut ref_file = File::open(REF_LOC).unwrap();
    let mut ref_string = String::new();
    let _ = ref_file.read_to_string(&mut ref_string);
    drop(ref_file);
    let ref_vec: Vec<&str> = ref_string.split('\n').collect();
    for ref_index in 0..(ref_vec.len()-1)/2 {
        let ref_name = &ref_vec[ref_index*2][1..];
        if !ref_name_vec.contains(&ref_name.to_string()) {
            ref_name_vec.push(ref_name.to_string());
            ref_struct_vec.push(
                Rc::new(Reference::new_raw(ref_name.to_string(), vec![Rc::new(Domain { 
                        cdd_acc: ref_name.split('|').last().unwrap().to_string(), 
                        super_acc: ref_name.split('|').last().unwrap().to_string()
                    })]))
            );
            let amr = ref_struct_vec.last().unwrap().get_classification();
            if db_amr_count.contains_key(&amr) {
                *db_amr_count.get_mut(&amr).unwrap() += 1;
            }
            else {
                db_amr_count.insert(amr, 1);
            }
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
    let mut amr_change_info: Vec<String> = vec![];

    // Make hashmap between reference name and reference struct to make matching constant time
    let ref_hashmap: StringHashMap<RefRc> = HashMap::from_iter(zip(ref_name_vec, ref_struct_vec));

    for sample_id in sample_id_vec {
    println!("{}", sample_id);
        for identity in vec![30, 50, 80] {
            println!("{}", identity);
            println!("LS");

            // Make distribution hashmap for diamond analysis
            let mut gene_amr_dia_hashmap: RcStringHashMap<usize> = HashMap::new();
            let mut cdd_amr_dia_hashmap: RcStringHashMap<usize> = HashMap::new();
            let mut super_amr_dia_hashmap: RcStringHashMap<usize> = HashMap::new();

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
                    if ls_query_struct_vec.len() > 0 {
                        // Count the last trio for diamond distr
                        ls_query_struct_vec.last().unwrap().borrow().update_amr_diamond_hashmaps(
                            &mut gene_amr_dia_hashmap, &mut cdd_amr_dia_hashmap, 
                            &mut super_amr_dia_hashmap, &mut domain_name_refs);
                    }
                    ls_query_struct_vec.push(RefCell::new(Query::new(&row, &ref_hashmap)));
                }
            }
            // Count the last trio for diamond distr
            ls_query_struct_vec.last().unwrap().borrow().update_amr_diamond_hashmaps(
                &mut gene_amr_dia_hashmap, &mut cdd_amr_dia_hashmap, 
                &mut super_amr_dia_hashmap, &mut domain_name_refs);

            // Make hashmap between query name and query struct to make matching constant time
            let ls_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(zip(ls_query_name_vec, ls_query_struct_vec));

            // Make hashmap for changes from one amr class to another
            let mut amr_change_hashmap: StringHashMap<usize> = HashMap::new();

            // Make distribution hashmap for deeparg analysis
            let mut gene_amr_dee_hashmap: RcStringHashMap<usize> = HashMap::new();
            let mut cdd_amr_dee_hashmap: RcStringHashMap<usize> = HashMap::new();
            let mut super_amr_dee_hashmap: RcStringHashMap<usize> = HashMap::new();

            // Make distribution hashmap for restricted diamond analysis
            let mut gene_amr_restr_dia_hashmap: RcStringHashMap<usize> = HashMap::new();
            let mut cdd_amr_restr_dia_hashmap: RcStringHashMap<usize> = HashMap::new();
            let mut super_amr_restr_dia_hashmap: RcStringHashMap<usize> = HashMap::new();

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

                // Count the switch
                let dee_amr = query_struct.get_top_deeparg_alignment().get_classification();
                let dia_amr = query_struct.get_top_diamond_alignment().get_classification();
                if dee_amr != dia_amr {
                    let dia_to_dee_amr = [dia_amr, dee_amr].join("\t");
                    if amr_change_hashmap.contains_key(&dia_to_dee_amr) {
                        *amr_change_hashmap.get_mut(&dia_to_dee_amr).unwrap() += 1;
                    }
                    else {
                        amr_change_hashmap.insert(dia_to_dee_amr, 1);
                    }
                }

                // Count the trio for deeparg and restricted diamond distr
                query_struct.update_amr_deeparg_hashmaps(
                    &mut gene_amr_dee_hashmap, &mut cdd_amr_dee_hashmap, 
                    &mut super_amr_dee_hashmap, &mut domain_name_refs);
                query_struct.update_amr_restricted_diamond_hashmaps(
                    &mut gene_amr_restr_dia_hashmap, &mut cdd_amr_restr_dia_hashmap, 
                    &mut super_amr_restr_dia_hashmap, &mut domain_name_refs);
            }

            // Save amr distribution info
            for (dia_to_dee_amr, count) in amr_change_hashmap.iter() {
                amr_change_info.push([sample_id.to_string(), format!("{}", identity), 
                    "LS".to_string(), dia_to_dee_amr.clone(), format!("{}", count)].join("\t"));
            }

            // Save trio and al. distribution info
            for query in ls_query_hashmap.values() {
                let borrowed = query.borrow();
                if !borrowed.is_deeparg_hit() {
                    continue;
                }
                if borrowed.are_diamond_and_deeparg_the_same() {
                    continue;
                }

                let (dia_gene_amr, dia_cdd_amr, dia_super_amr) =
                    borrowed.get_top_diamond_alignment_domain_identifiers();
                let (dee_gene_amr, dee_cdd_amr, dee_super_amr) = 
                    borrowed.get_top_deeparg_hit_domain_identifiers().unwrap();
                if dia_gene_amr == dee_gene_amr {
                    continue;
                }
                let dee_gene_amr_info = 
                    get_domain_info(dee_gene_amr.to_string(), &gene_amr_dee_hashmap,
                        &gene_amr_restr_dia_hashmap, &gene_amr_dia_hashmap, 
                        &db_gene_amr_hashmap, &mut domain_name_refs);
                let dee_cdd_amr_info = 
                    get_domain_info(dee_cdd_amr.to_string(), &cdd_amr_dee_hashmap,
                        &cdd_amr_restr_dia_hashmap, &cdd_amr_dia_hashmap, 
                        &db_cdd_amr_hashmap, &mut domain_name_refs);
                let dee_super_amr_info = 
                    get_domain_info(dee_super_amr.to_string(), &super_amr_dee_hashmap,
                        &super_amr_restr_dia_hashmap, &super_amr_dia_hashmap, 
                        &db_super_amr_hashmap, &mut domain_name_refs);
                let dia_gene_amr_info = 
                    get_domain_info(dia_gene_amr.to_string(), &gene_amr_dee_hashmap,
                        &gene_amr_restr_dia_hashmap, &gene_amr_dia_hashmap, 
                        &db_gene_amr_hashmap, &mut domain_name_refs);
                let dia_cdd_amr_info = 
                    get_domain_info(dia_cdd_amr.to_string(), &cdd_amr_dee_hashmap,
                        &cdd_amr_restr_dia_hashmap, &cdd_amr_dia_hashmap, 
                        &db_cdd_amr_hashmap, &mut domain_name_refs);
                let dia_super_amr_info = 
                    get_domain_info(dia_super_amr.to_string(), &super_amr_dee_hashmap,
                        &super_amr_restr_dia_hashmap, &super_amr_dia_hashmap, 
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

            println!("SS");
            // Clear distribution hashmap for short sequence diamond analysis
            gene_amr_dia_hashmap.clear();
            cdd_amr_dia_hashmap.clear();
            super_amr_dia_hashmap.clear();

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
                    if ss_query_struct_vec.len() > 0 {
                        // Count the last trio for diamond distr
                        ss_query_struct_vec.last().unwrap().borrow().update_amr_diamond_hashmaps(
                            &mut gene_amr_dia_hashmap, &mut cdd_amr_dia_hashmap, 
                            &mut super_amr_dia_hashmap, &mut domain_name_refs);
                    }
                    ss_query_struct_vec.push(RefCell::new(Query::new(&row, &ref_hashmap)));
                    
                }
            }
            // Count the last trio for diamond distr
            ss_query_struct_vec.last().unwrap().borrow().update_amr_diamond_hashmaps(
                    &mut gene_amr_dia_hashmap, &mut cdd_amr_dia_hashmap, 
                    &mut super_amr_dia_hashmap, &mut domain_name_refs);

            // Make hashmap between query name and query struct to make matching constant time
            let ss_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(
                zip(ss_query_name_vec, ss_query_struct_vec));

            // Clear hashmap for changes from one amr class to another for short sequence 
            amr_change_hashmap.clear();

            // Clear distribution hashmap for short sequence deeparg analysis
            gene_amr_dee_hashmap.clear();
            cdd_amr_dee_hashmap.clear();
            super_amr_dee_hashmap.clear();

            // Clear distribution hashmap for short sequence restricted diamond analysis
            gene_amr_restr_dia_hashmap.clear();
            cdd_amr_restr_dia_hashmap.clear();
            super_amr_restr_dia_hashmap.clear();

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

                // Count the switch
                let dee_amr = query_struct.get_top_deeparg_alignment().get_classification();
                let dia_amr = query_struct.get_top_diamond_alignment().get_classification();
                if dee_amr != dia_amr {
                    let dia_to_dee_amr = [dia_amr, dee_amr].join("\t");
                    if amr_change_hashmap.contains_key(&dia_to_dee_amr) {
                        *amr_change_hashmap.get_mut(&dia_to_dee_amr).unwrap() += 1;
                    }
                    else {
                        amr_change_hashmap.insert(dia_to_dee_amr, 1);
                    }
                }

                // Count the trio for deeparg and restricted diamond distr
                query_struct.update_amr_deeparg_hashmaps(
                    &mut gene_amr_dee_hashmap, &mut cdd_amr_dee_hashmap, 
                    &mut super_amr_dee_hashmap, &mut domain_name_refs);
                query_struct.update_amr_restricted_diamond_hashmaps(
                    &mut gene_amr_restr_dia_hashmap, &mut cdd_amr_restr_dia_hashmap, 
                    &mut super_amr_restr_dia_hashmap, &mut domain_name_refs);
            }

            // Save amr distribution info
            for (dia_to_dee_amr, count) in amr_change_hashmap.iter() {
                amr_change_info.push([sample_id.to_string(), format!("{}", identity), 
                    "SS".to_string(), dia_to_dee_amr.clone(), format!("{}", count)].join("\t"));
            }

            // Save trio and al. distribution info
            for query in ss_query_hashmap.values() {
                let borrowed = query.borrow();
                if !borrowed.is_deeparg_hit() {
                    continue;
                }
                if borrowed.are_diamond_and_deeparg_the_same() {
                    continue;
                }

                let (dia_gene_amr, dia_cdd_amr, dia_super_amr) = 
                    borrowed.get_top_diamond_alignment_domain_identifiers();
                let (dee_gene_amr, dee_cdd_amr, dee_super_amr) = 
                    borrowed.get_top_deeparg_hit_domain_identifiers().unwrap();
                if dia_gene_amr == dee_gene_amr {
                    continue;
                }
                let dee_gene_amr_info = 
                    get_domain_info(dee_gene_amr.to_string(), &gene_amr_dee_hashmap,
                        &gene_amr_restr_dia_hashmap, &gene_amr_dia_hashmap, 
                        &db_gene_amr_hashmap, &mut domain_name_refs);
                let dee_cdd_amr_info = 
                    get_domain_info(dee_cdd_amr.to_string(), &cdd_amr_dee_hashmap,
                        &cdd_amr_restr_dia_hashmap, &cdd_amr_dia_hashmap, 
                        &db_cdd_amr_hashmap, &mut domain_name_refs);
                let dee_super_amr_info = 
                    get_domain_info(dee_super_amr.to_string(), &super_amr_dee_hashmap,
                        &super_amr_restr_dia_hashmap, &super_amr_dia_hashmap, 
                        &db_super_amr_hashmap, &mut domain_name_refs);
                let dia_gene_amr_info = 
                    get_domain_info(dia_gene_amr.to_string(), &gene_amr_dee_hashmap,
                        &gene_amr_restr_dia_hashmap, &gene_amr_dia_hashmap, 
                        &db_gene_amr_hashmap, &mut domain_name_refs);
                let dia_cdd_amr_info = 
                    get_domain_info(dia_cdd_amr.to_string(), &cdd_amr_dee_hashmap,
                        &cdd_amr_restr_dia_hashmap, &cdd_amr_dia_hashmap, 
                        &db_cdd_amr_hashmap, &mut domain_name_refs);
                let dia_super_amr_info = 
                    get_domain_info(dia_super_amr.to_string(), &super_amr_dee_hashmap,
                        &super_amr_restr_dia_hashmap, &super_amr_dia_hashmap, 
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

    // Create output files
    let mut comparison_output = File::create(COMPARISON_OUTPUT).unwrap();
    let _ = comparison_output.write(["Sample", "Alignment Identity", "Model", 
        "Diamond alignment cdd|gene|amr key", "Reference distribution",
        "Diamond distribution", "Diamond restricted distribution", 
        "Deeparg distribution", "Diamond alignment cdd|amr key",
        "Reference distribution", "Diamond distribution", 
        "Diamond restricted distribution", "Deeparg distribution", 
        "Diamond alignment super|amr key", "Reference distribution",
        "Diamond distribution", "Diamond restricted distribution", 
        "Deeparg distribution", "DeepARG hit cdd|gene|amr key",
        "Reference distribution", "Diamond distribution", 
        "Diamond restricted distribution", "Deeparg distribution",
        "DeepARG hit cdd|amr key", "Reference distribution",
        "Diamond distribution", "Diamond restricted distribution", 
        "Deeparg distribution", "DeepARG hit super|amr key",
        "Reference distribution", "Diamond distribution", 
        "Diamond restricted distribution", "Deeparg distribution", 
        "Diamond alignment cdd|gene|amr/DeepARG hit cdd|gene|amr deeparg count\n"
        ].join("\t").as_bytes());
    write_hashmap_to_file(&mut comparison_output, deeparg_not_equal_diamond);
    drop(comparison_output);

    let mut amr_switch_output = File::create(AMR_SWITCH_DISTRIBUTION).unwrap();
    let _ = amr_switch_output.write(b"AMR classification\tCount\n");
    for (amr_class, count) in db_amr_count.iter() {
        let _ = amr_switch_output.write(format!("{}\t{}\n", amr_class, count).as_bytes());
    }
    let _ = amr_switch_output.write(b"Sample\tAlignment Identity\tModel\tDiamond AMR\tDeepARG AMR\tCount\n");
    for row in amr_change_info.iter() {
        let _ = amr_switch_output.write(format!("{}\n", row).as_bytes());
    }
    drop(amr_switch_output);
}

fn get_domain_info(dom: String, 
        deeparg_hashmap: &RcStringHashMap<usize>, restricted_diamond_hashmap: &RcStringHashMap<usize>, 
        diamond_hashmap: &RcStringHashMap<usize>, hashmap_db: &RcStringHashMap<usize>, 
        domain_name_refs: &mut StringHashMap<Rc<String>>) -> String {
    let dom_refs = domain_name_refs.get(&dom).unwrap();
    if restricted_diamond_hashmap.get(dom_refs).is_none() && deeparg_hashmap.get(dom_refs).is_none() {
        panic!("Messed up with {dom}")
    }
    format!("{}\t{}\t{}\t{}\t{}", dom.clone(), *hashmap_db.get(dom_refs).unwrap(), 
        *diamond_hashmap.get(dom_refs).unwrap_or(&0), *restricted_diamond_hashmap.get(dom_refs).unwrap_or(&0),
        *deeparg_hashmap.get(dom_refs).unwrap_or(&0))
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

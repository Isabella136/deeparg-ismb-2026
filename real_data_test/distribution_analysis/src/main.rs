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

static CDD_DIR: &str = "../../CDD_features/";
static REF_LOC: &str = "../../data/database/v2/features.fasta";
static DEEPARG_HIT_FILE: &str = "X.mapping.ARG";
static ALIGNMENT_FILE: &str = "X.align.daa.tsv";
static SAMPLE_ID_FILE: &str = "../real_samples.txt";
static PATH_TO_SS: &str = "deeparg_results";
static PATH_TO_LS: &str = "spades/deeparg_results";

static DISTRIBUTION_OUTPUT: &str = "_distributions.txt";

fn main() {
    // Read content of real_samples.txt to find biosample ID
    let mut sample_id_file = File::open(SAMPLE_ID_FILE).unwrap();
    let mut sample_id_string = String::new();
    let _ = sample_id_file.read_to_string(&mut sample_id_string);
    drop(sample_id_file);
    let sample_id_vec: Vec<&str> = sample_id_string.split('\n').collect();

    // Make distribution hashmap for entire database
    let mut db_gene_amr_hashmap: StringHashMap<f64> = HashMap::new();
    let mut db_cdd_amr_hashmap: StringHashMap<f64> = HashMap::new();
    let mut db_super_amr_hashmap: StringHashMap<f64> = HashMap::new();

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
                        &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, &mut db_super_amr_hashmap);
                }
                ref_name_vec.push(name.to_string());
                ref_struct_vec.push(Rc::new(Reference::new(row)));
            }
        }
    }

    ref_struct_vec.last().unwrap().update_amr_hashmaps(
        &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, &mut db_super_amr_hashmap);
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
                &mut db_gene_amr_hashmap, &mut db_cdd_amr_hashmap, &mut db_super_amr_hashmap);
        }
    }

    // Create tables that includes db distribution followed by [30SS, 30LS, 50SS, 50LS, 80SS, 80LS]*10 distributions
    let mut gene_amr_hashtable_a1: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_gene_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));
    let mut cdd_amr_hashtable_a1: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_cdd_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));
    let mut super_amr_hashtable_a1: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_super_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));

    let mut gene_amr_hashtable_a2: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_gene_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));
    let mut cdd_amr_hashtable_a2: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_cdd_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));
    let mut super_amr_hashtable_a2: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_super_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));

    let mut gene_amr_hashtable_a3: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_gene_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));
    let mut cdd_amr_hashtable_a3: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_cdd_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));
    let mut super_amr_hashtable_a3: StringHashMap<Vec<f64>> = HashMap::from_iter(
        db_super_amr_hashmap.iter().map(|(x, y)| (x.clone(), vec![*y])));

    // Make hashmap between reference name and reference struct to make matching constant time
    let ref_hashmap: StringHashMap<RefRc> = HashMap::from_iter(zip(ref_name_vec, ref_struct_vec));

    for sample_id in sample_id_vec {
        for identity in vec![30, 50, 80] {
            // Make distribution hashmap for analysis 3
            let mut gene_amr_hashmap3: StringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap3: StringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap3: StringHashMap<f64> = HashMap::new();

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
                        &mut gene_amr_hashmap3, &mut cdd_amr_hashmap3, &mut super_amr_hashmap3);
                }
            }
            // Make hashmap between query name and query struct to make matching constant time
            let ss_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(
                zip(ss_query_name_vec, ss_query_struct_vec));

            // Make distribution hashmap for analysis 1
            let mut gene_amr_hashmap1: StringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap1: StringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap1: StringHashMap<f64> = HashMap::new();

            // Make distribution hashmap for analysis 2
            let mut gene_amr_hashmap2: StringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap2: StringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap2: StringHashMap<f64> = HashMap::new();

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
                    &mut gene_amr_hashmap1, &mut cdd_amr_hashmap1, &mut super_amr_hashmap1);
                query_struct.update_amr_hashmaps_v2(
                    &mut gene_amr_hashmap2, &mut cdd_amr_hashmap2, &mut super_amr_hashmap2);
            }

            gene_amr_hashtable_a1 = return_extended_hashtable(&mut gene_amr_hashtable_a1, 
                gene_amr_hashmap1);
            cdd_amr_hashtable_a1 = return_extended_hashtable(&mut cdd_amr_hashtable_a1, 
                cdd_amr_hashmap1);
            super_amr_hashtable_a1 = return_extended_hashtable(&mut super_amr_hashtable_a1, 
                super_amr_hashmap1);

            gene_amr_hashtable_a2 = return_extended_hashtable(&mut gene_amr_hashtable_a2, 
                gene_amr_hashmap2);
            cdd_amr_hashtable_a2 = return_extended_hashtable(&mut cdd_amr_hashtable_a2, 
                cdd_amr_hashmap2);
            super_amr_hashtable_a2 = return_extended_hashtable(&mut super_amr_hashtable_a2, 
                super_amr_hashmap2);

            gene_amr_hashtable_a3 = return_extended_hashtable(&mut gene_amr_hashtable_a3, 
                gene_amr_hashmap3);
            cdd_amr_hashtable_a3 = return_extended_hashtable(&mut cdd_amr_hashtable_a3, 
                cdd_amr_hashmap3);
            super_amr_hashtable_a3 = return_extended_hashtable(&mut super_amr_hashtable_a3, 
                super_amr_hashmap3);

            // Make distribution hashmap for analysis 3
            let mut gene_amr_hashmap3: StringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap3: StringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap3: StringHashMap<f64> = HashMap::new();

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
                        &mut gene_amr_hashmap3, &mut cdd_amr_hashmap3, &mut super_amr_hashmap3);
                }
            }
            // Make hashmap between query name and query struct to make matching constant time
            let ls_query_hashmap: StringHashMap<QueRef> = HashMap::from_iter(zip(ls_query_name_vec, ls_query_struct_vec));

            // Make distribution hashmap for analysis 1
            let mut gene_amr_hashmap1: StringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap1: StringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap1: StringHashMap<f64> = HashMap::new();

            // Make distribution hashmap for analysis 2
            let mut gene_amr_hashmap2: StringHashMap<f64> = HashMap::new();
            let mut cdd_amr_hashmap2: StringHashMap<f64> = HashMap::new();
            let mut super_amr_hashmap2: StringHashMap<f64> = HashMap::new();

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
                    &mut gene_amr_hashmap1, &mut cdd_amr_hashmap1, &mut super_amr_hashmap1);
                query_struct.update_amr_hashmaps_v2(
                    &mut gene_amr_hashmap2, &mut cdd_amr_hashmap2, &mut super_amr_hashmap2);
            }

            gene_amr_hashtable_a1 = return_extended_hashtable(&mut gene_amr_hashtable_a1, 
                gene_amr_hashmap1);
            cdd_amr_hashtable_a1 = return_extended_hashtable(&mut cdd_amr_hashtable_a1, 
                cdd_amr_hashmap1);
            super_amr_hashtable_a1 = return_extended_hashtable(&mut super_amr_hashtable_a1, 
                super_amr_hashmap1);

            gene_amr_hashtable_a2 = return_extended_hashtable(&mut gene_amr_hashtable_a2, 
                gene_amr_hashmap2);
            cdd_amr_hashtable_a2 = return_extended_hashtable(&mut cdd_amr_hashtable_a2, 
                cdd_amr_hashmap2);
            super_amr_hashtable_a2 = return_extended_hashtable(&mut super_amr_hashtable_a2, 
                super_amr_hashmap2);

            gene_amr_hashtable_a3 = return_extended_hashtable(&mut gene_amr_hashtable_a3, 
                gene_amr_hashmap3);
            cdd_amr_hashtable_a3 = return_extended_hashtable(&mut cdd_amr_hashtable_a3, 
                cdd_amr_hashmap3);
            super_amr_hashtable_a3 = return_extended_hashtable(&mut super_amr_hashtable_a3, 
                super_amr_hashmap3);
        }
    }

    // Create output files
    let mut gene_amr_distr_a1 = File::create(format!("analysis_1_gene_amr{DISTRIBUTION_OUTPUT}")).unwrap();
    let mut cdd_amr_distr_a1 = File::create(format!("analysis_1_cdd_amr{DISTRIBUTION_OUTPUT}")).unwrap();
    let mut super_amr_distr_a1 = File::create(format!("analysis_1_super_amr{DISTRIBUTION_OUTPUT}")).unwrap();

    let mut gene_amr_distr_a2 = File::create(format!("analysis_2_gene_amr{DISTRIBUTION_OUTPUT}")).unwrap();
    let mut cdd_amr_distr_a2 = File::create(format!("analysis_2_cdd_amr{DISTRIBUTION_OUTPUT}")).unwrap();
    let mut super_amr_distr_a2 = File::create(format!("analysis_2_super_amr{DISTRIBUTION_OUTPUT}")).unwrap();

    let mut gene_amr_distr_a3 = File::create(format!("analysis_3_gene_amr{DISTRIBUTION_OUTPUT}")).unwrap();
    let mut cdd_amr_distr_a3 = File::create(format!("analysis_3_cdd_amr{DISTRIBUTION_OUTPUT}")).unwrap();
    let mut super_amr_distr_a3 = File::create(format!("analysis_3_super_amr{DISTRIBUTION_OUTPUT}")).unwrap();



    write_hashtable_to_file(&mut gene_amr_distr_a1, gene_amr_hashtable_a1);
    write_hashtable_to_file(&mut cdd_amr_distr_a1, cdd_amr_hashtable_a1);
    write_hashtable_to_file(&mut super_amr_distr_a1, super_amr_hashtable_a1);

    write_hashtable_to_file(&mut gene_amr_distr_a2, gene_amr_hashtable_a2);
    write_hashtable_to_file(&mut cdd_amr_distr_a2, cdd_amr_hashtable_a2);
    write_hashtable_to_file(&mut super_amr_distr_a2, super_amr_hashtable_a2);

    write_hashtable_to_file(&mut gene_amr_distr_a3, gene_amr_hashtable_a3);
    write_hashtable_to_file(&mut cdd_amr_distr_a3, cdd_amr_hashtable_a3);
    write_hashtable_to_file(&mut super_amr_distr_a3, super_amr_hashtable_a3);
}

fn write_hashtable_to_file(file: &mut File, table: StringHashMap<Vec<f64>>) {
    for row in table.iter() {
        let _ = file.write(row.0.as_bytes());
        for val in row.1 {
            let _ = file.write(format!("\t{}", *val).as_bytes());
        }
        let _ = file.write(b"\n");
    }
    

}

fn return_extended_hashtable(hashtable: &mut StringHashMap<Vec<f64>>, 
        hashmap: StringHashMap<f64>) -> StringHashMap<Vec<f64>>{
    HashMap::from_iter(hashtable.iter_mut().map(|(x, y)| {
        if hashmap.contains_key(x) {y.push(*hashmap.get(x).unwrap());}
        else {y.push(0.0);}
        (x.clone(), y.clone())
    }))
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

    fn update_amr_hashmaps(&self, gene_map: &mut StringHashMap<f64>, 
            cdd_map: &mut StringHashMap<f64>, super_map: &mut StringHashMap<f64>){
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
                domain.cdd_acc, self.get_gene_name(), self.get_classification());
            if gene_map.contains_key(&gene_key) {
                *gene_map.get_mut(&gene_key).unwrap() += weight
                
            }
            else {
                gene_map.insert(gene_key, weight);
            }
            let cdd_key = format!("{}|{}", 
                domain.cdd_acc, self.get_classification());
            if cdd_map.contains_key(&cdd_key) {
                *cdd_map.get_mut(&cdd_key).unwrap() += weight;
            }
            else {
                cdd_map.insert(cdd_key, weight);
            }
            let super_key = format!("{}|{}", 
                domain.super_acc, self.get_classification());
            if super_map.contains_key(&super_key) {
                *super_map.get_mut(&super_key).unwrap() += weight;
            }
            else {
                super_map.insert(super_key, weight);
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
        let mut domains = vec![];
        for domain in matching_reference.domains.iter() {
            if (domain.start >= start && domain.start < end) || (domain.end <= end && domain.end > start) {
                domains.push(domain.clone())
            }
        }
        Self { matching_reference, start, end, domains}
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
    fn add_alignment(&mut self, row: &Vec<&str>, ref_hashmap: &StringHashMap<RefRc>) {
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
    fn update_amr_hashmaps_v1(&self, gene_map: &mut StringHashMap<f64>, 
            cdd_map: &mut StringHashMap<f64>, super_map: &mut StringHashMap<f64>) {
        if self.top_deeparg_hit.is_some() {
            self.alignments[self.top_deeparg_hit.unwrap()]
                .update_amr_hashmaps(gene_map, cdd_map, super_map);
        }
    }
    fn update_amr_hashmaps_v2(&self, gene_map: &mut StringHashMap<f64>, 
            cdd_map: &mut StringHashMap<f64>, super_map: &mut StringHashMap<f64>) {
        if self.top_deeparg_hit.is_some() {
            self.alignments[self.top_diamond_alignment]
                .update_amr_hashmaps(gene_map, cdd_map, super_map);
        }
    }
    fn update_amr_hashmaps_v3(&self, gene_map: &mut StringHashMap<f64>, 
            cdd_map: &mut StringHashMap<f64>, super_map: &mut StringHashMap<f64>) {
        self.alignments[self.top_diamond_alignment]
            .update_amr_hashmaps(gene_map, cdd_map, super_map);
    }
}
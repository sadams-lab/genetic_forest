use std::fs::File;
use std::str;
use sprs::{CsMat, Shape, TriMat};
use rand::{Rng, thread_rng};
use rand::seq::SliceRandom;

pub struct GenoMatrix {
    pub ids: Vec<String>,
    pub phenotypes: Vec<u64>,
    pub n_subjects: f64,
    pub n_genotypes: f64,
    pub n_genotypes_minus_mask: f64,
    var_mask: Vec<usize>,
    pheno_weight: f64, // weight to use for selecting phenotype (e.g. 0.5 would be balanced...)
    genotypes: CsMat<u8>
}

pub struct GenoMatrixSlice {
    pub subj_ids: Vec<usize>,
    pub genotype_ids: Vec<usize>,
}

impl GenoMatrix {
    
    pub fn new(path: &str, sep: &str) -> Self {
        let mut row_ids: Vec<String> = Vec::new();
        let mut phenotypes: Vec<u64> = Vec::new();
        let mat_size: Shape = get_mat_size(path, sep);
        let mut geno_mat = TriMat::new(mat_size);
        let mut rdr = make_reader(path, sep);
        let mut rownum = 0;
        for result in rdr.records() {
            let mut colnum = 2;
            let record = result.unwrap();
            row_ids.push(record[0].to_string());
            phenotypes.push(record[1].parse::<u64>().unwrap());
            loop {
                if colnum == mat_size.1 + 2 {
                    break
                }
                geno_mat.add_triplet(rownum, colnum - 2, record[colnum].parse::<u8>().unwrap());
                colnum += 1;
            }
            rownum += 1
        }
        let pheno_weight: f64 = phenotypes.iter().sum::<u64>() as f64 / mat_size.0 as f64;
        GenoMatrix{
            ids: row_ids, 
            phenotypes: phenotypes,
            n_subjects: mat_size.0 as f64,
            n_genotypes: mat_size.1 as f64,
            n_genotypes_minus_mask: mat_size.1 as f64,
            var_mask: Vec::new(),
            genotypes: geno_mat.to_csr(),
            pheno_weight: pheno_weight
        }
    }

    pub fn make_slice(&self, n_vars: &f64, n_subj: &f64) -> GenoMatrixSlice {
/*         Note that these are all sampling without replacement
        We don't implement sampling with replacement for this
        Sampling with replacement improves predictive ability of the model
        Which we do not care about 
        
        Subject selection is auto-weighted for phenotype. 
        Might consider adding an option to disable this or accept 
        external weights (or just specify target ratio)*/
        let select_var_prob: f64 = n_vars / self.n_genotypes_minus_mask;
        let select_subj_prob: f64 = n_subj / self.n_subjects;
        let prob_0 = self.pheno_weight * select_subj_prob;
        let prob_1 = (1. - self.pheno_weight) * select_subj_prob;
        let mut subjs: Vec<usize> = Vec::new();
        let mut g_ids: Vec<usize> = Vec::new();
        let mut rng = thread_rng();
        for s in 0..self.n_subjects as usize {
            match self.phenotypes[s] {
                0 => {
                    if rng.gen_bool(prob_0) {
                        subjs.push(s);
                    }
                },
                1 => {
                    if rng.gen_bool(prob_1) {
                        subjs.push(s);
                    }
                },
                _ => ()
            }
        };
        for g in 0..self.n_genotypes as usize {
            if self.var_mask.iter().any(|i| i == &g) {
                (); // skip variants that are masked
            } 
            else {
                if rng.gen_bool(select_var_prob) {
                    g_ids.push(g);
                };
            };
        };
        GenoMatrixSlice {
            subj_ids: subjs,
            genotype_ids: g_ids
        }
    }
    
    pub fn get_slice_data(&self, gm: &GenoMatrixSlice) -> (Vec<&u64>, Vec<&u64>, Vec<Vec<&u8>>) {
        let mut rng = thread_rng();
        let mut g_vec: Vec<Vec<&u8>> = Vec::new();
        let mut p_vec: Vec<&u64> = Vec::new();
        for s in &gm.subj_ids {
            p_vec.push(&self.phenotypes[*s]);
        };
        for g in &gm.genotype_ids {
            let mut gv: Vec<&u8> = Vec::new();
            for s in &gm.subj_ids {
                gv.push(self.genotypes.get(*s, *g).unwrap());
            }
            g_vec.push(gv);
        };
        let mut pheno2: Vec<&u64> = p_vec.to_vec();
        pheno2.shuffle(&mut rng);
        return (p_vec, pheno2, g_vec)
    }

    pub fn mask_vars(&mut self, variants: Vec<usize>) {
        // add variants to mask from analysis
        // also adjust the number of variants
        self.n_genotypes_minus_mask -= variants.len() as f64;
        for v in variants {
            self.var_mask.push(v);
        }
    }
}

fn make_reader(path: &str, sep: &str) -> csv::Reader<File> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(why) => panic!("Something happened with the file: {}", why),
    };
    return csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(sep.as_bytes()[0])
        .comment(Some(b'#'))
        .from_reader(file);
}

fn get_mat_size(path: &str, sep: &str) -> Shape {
    let mut ncols: usize = 0;
    let mut nrows: usize = 0;
    let mut rdr = make_reader(path, sep);
    for result in rdr.byte_records() {
            let record = result.unwrap();
            if ncols == 0 {
                ncols = &record.len() - 2;
            };
            nrows += 1;
    }
    (nrows, ncols)
}

use sprs::{CsMat, Shape, TriMat};
use rand::{Rng, thread_rng};
use rand::seq::SliceRandom;
use std::fs::File;

pub struct GenoMatrix {
    pub ids: Vec<String>,
    pub phenotypes: Vec<f64>,
    pub n_subjects: f64,
    pub n_genotypes: f64,
    pheno_weight: f64, // weight to use for selecting phenotype (e.g. 0.5 would be balanced...)
    genotypes: CsMat<u8>,
    genotype_indices: Vec<usize>
}

pub struct GenoMatrixSlice {
    pub subj_ids: Vec<usize>,
    pub genotype_ids: Vec<usize>,
}

impl GenoMatrix {
    
    pub fn new(rdr: &mut csv::Reader<File>, mat_size: Shape, continuous_outcome: &bool) -> Self {
        let pheno_weight: f64;
        let mut row_ids: Vec<String> = Vec::new();
        let mut phenotypes: Vec<f64> = Vec::new();
        let mut geno_mat = TriMat::new(mat_size);
        let mut rownum = 0;
        for result in rdr.records() {
            let mut colnum = 2;
            let record = result.unwrap();
            row_ids.push(record[0].to_string());
            phenotypes.push(record[1].parse::<f64>().unwrap());
            loop {
                if colnum == mat_size.1 + 2 {
                    break
                }
                geno_mat.add_triplet(rownum, colnum - 2, record[colnum].parse::<u8>().unwrap());
                colnum += 1;
            }
            rownum += 1
        }
        if *continuous_outcome {
            pheno_weight = -1.
        } else {
            pheno_weight = phenotypes.iter().sum::<f64>() as f64 / mat_size.0 as f64;
        }
        GenoMatrix{
            ids: row_ids, 
            phenotypes: phenotypes,
            n_subjects: mat_size.0 as f64,
            n_genotypes: mat_size.1 as f64,
            genotypes: geno_mat.to_csr(),
            pheno_weight: pheno_weight,
            genotype_indices: (0..mat_size.1).collect()
        }
    }

    pub fn make_slice(&self, var_frac: f64, subj_frac: f64) -> GenoMatrixSlice {
    //  Note that these are all sampling without replacement
    //  We don't implement sampling with replacement for this
    //  Sampling with replacement improves predictive ability of the model
    //  Which we do not care about.
    //  Subject selection is auto-weighted for phenotype. 
    //  Might consider adding an option to disable this or accept 
    //  external weights (or just specify target ratio)

        let prob_0 = self.pheno_weight * subj_frac;
        let prob_1 = (1. - self.pheno_weight) * subj_frac;
        let mut subjs: Vec<usize> = Vec::new();
        let mut g_ids: Vec<usize> = Vec::new();
        let mut rng = thread_rng();
        for s in 0..self.n_subjects as usize {
            match self.phenotypes[s] {
                _x if self.pheno_weight == -1.0 => { // happens if pheno is continuous, no weighting applied
                    if rng.gen_bool(subj_frac) {
                        subjs.push(s);
                    }
                },
                _x if _x == 0.0 => {
                    if rng.gen_bool(prob_0) {
                        subjs.push(s);
                    }
                },
                _x if _x == 1.0 => {
                    if rng.gen_bool(prob_1) {
                        subjs.push(s);
                    }
                },
                _ => ()
            }
        };
        for g in &self.genotype_indices {
            if rng.gen_bool(var_frac) {
                g_ids.push(*g);
            };
        };
        GenoMatrixSlice {
            subj_ids: subjs,
            genotype_ids: g_ids
        }
    }
    
    pub fn get_slice_data(&self, gm: &GenoMatrixSlice) -> (Vec<&f64>, Vec<&f64>, Vec<Vec<&u8>>) {
        let mut rng = thread_rng();
        let mut g_vec: Vec<Vec<&u8>> = Vec::new();
        let mut p_vec: Vec<&f64> = Vec::new();
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
        let mut pheno2: Vec<&f64> = p_vec.to_vec();
        pheno2.shuffle(&mut rng);
        return (p_vec, pheno2, g_vec)
    }

    pub fn set_genotype_indices(&mut self, variants: Vec<usize>) {
        self.genotype_indices = variants;
    }
}

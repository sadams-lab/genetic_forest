use std::fs::File;
use std::str;
use sprs::{CsMat, Shape, TriMat};
use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;

pub struct GenoMatrix {
    pub ids: Vec<String>,
    pub phenotypes: Vec<u8>,
    pub n_subjects: f64,
    pub n_genotypes: f64,
    genotypes: CsMat<u8>
}

pub struct GenoMatrixSlice {
    pub subj_ids: Vec<usize>,
    pub genotype_ids: Vec<usize>,
}

impl GenoMatrix {
    
    pub fn new(path: &str) -> Self {
        let mut row_ids: Vec<String> = Vec::new();
        let mut phenotypes: Vec<u8> = Vec::new();
        let mat_size: Shape = get_mat_size(path);
        let mut geno_mat = TriMat::new(mat_size);
        let mut rdr = make_reader(path);
        let mut rownum = 0;
        for result in rdr.records() {
            let mut colnum = 2;
            let record = result.unwrap();
            row_ids.push(record[0].to_string());
            phenotypes.push(record[1].parse::<u8>().unwrap());
            loop {
                if colnum == mat_size.1 + 2 {
                    break
                }
                geno_mat.add_triplet(rownum, colnum - 2, record[colnum].parse::<u8>().unwrap());
                colnum += 1;
            }
            rownum += 1
        }
        GenoMatrix{
            ids: row_ids, 
            phenotypes: phenotypes,
            n_subjects: mat_size.0 as f64,
            n_genotypes: mat_size.1 as f64,
            genotypes: geno_mat.to_csr(),
        }
    }

    pub fn make_slice(&self, n_vars: &f64, n_subj: &f64) -> GenoMatrixSlice {
        let select_var_prob: f64 = n_vars / self.n_genotypes;
        let select_subj_prob: f64 = n_subj / self.n_subjects;
        let mut subjs: Vec<usize> = Vec::new();
        let mut g_ids: Vec<usize> = Vec::new();
        let mut rng = thread_rng();
        for s in 0..self.n_subjects as usize {
            if rng.gen_bool(select_subj_prob) {
                subjs.push(s);

            };
        };
        for g in 0..self.n_genotypes as usize {
            if rng.gen_bool(select_var_prob) {
                g_ids.push(g);
            };
        };
        GenoMatrixSlice {
            subj_ids: subjs,
            genotype_ids: g_ids
        }
    }
    
    pub fn get_slice_data(&self, gm: &GenoMatrixSlice) -> (Vec<&u8>, Vec<&u8>, Vec<Vec<&u8>>) {
        let mut g_vec: Vec<Vec<&u8>> = Vec::new();
        let mut p_vec: Vec<&u8> = Vec::new();
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
        let mut pheno2: Vec<&u8> = p_vec.to_vec();
        pheno2.shuffle(&mut thread_rng());
        return (p_vec, pheno2, g_vec)
    }
}

fn make_reader(path: &str) -> csv::Reader<File> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(why) => panic!("Something happened with the file: {}", why),
    };
    return csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .from_reader(file);
}

fn get_mat_size(path: &str) -> Shape {
    let mut ncols: usize = 0;
    let mut nrows: usize = 0;
    let mut rdr = make_reader(path);
    for result in rdr.byte_records() {
            let record = result.unwrap();
            if ncols == 0 {
                ncols = &record.len() - 2;
            };
            nrows += 1;
    }
    (nrows, ncols)
}

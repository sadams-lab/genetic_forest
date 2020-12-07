use crate::tree;
use crate::matrix;
use rayon::prelude::*;

pub struct Forest {
    pub trees: Vec<tree::Node>
}

impl Forest {
    pub fn grow(gm: matrix::GenoMatrix, n_tree: i32, mtry: f64, min_node_size: i32, subj_fraction: f64, threads: usize) -> Self {
        let var_sample_size: f64 = mtry * gm.n_genotypes;
        let subj_sample_size: f64 = subj_fraction * gm.n_subjects;
        //let mut trees: Vec<tree::Node> = Vec::new();
        //let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let t: Vec<tree::Node> = (0..n_tree).into_par_iter()
        .map(|i| -> tree::Node {
            make_tree(&gm, &var_sample_size, &subj_sample_size, &min_node_size)
        }).collect();
        Forest {trees: t}
    }
}

fn make_tree(gm: &matrix::GenoMatrix, n_vars: &f64, n_subj: &f64, min_node_size: &i32) -> tree::Node {
    let sample = gm.make_slice(n_vars, n_subj);
    let data = gm.get_slice_data(&sample);
    tree::Node::grow(data, *min_node_size, sample)
}
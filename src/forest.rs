use crate::tree;
use crate::matrix;
use rayon::prelude::*;

use std::collections::HashMap;


pub struct Forest {
    pub trees: Vec<tree::Node>
}

impl Forest {
    pub fn grow(gm: matrix::GenoMatrix, n_tree: i32, mtry: f64, min_node_size: i32, subj_fraction: f64, threads: usize, random_seed: u64) -> Self {
        let var_sample_size: f64 = mtry * gm.n_genotypes;
        let subj_sample_size: f64 = subj_fraction * gm.n_subjects;
        //let mut trees: Vec<tree::Node> = Vec::new();
        //let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let t: Vec<tree::Node> = (0..n_tree).into_par_iter()
        .map(|i| -> tree::Node {
            make_tree(&gm, &var_sample_size, &subj_sample_size, &min_node_size, &random_seed)
        }).collect();
        Forest {trees: t}
    }
    
    pub fn var_importance(&self) {
        let mut tree_imps: HashMap<usize, Vec<f32>> = HashMap::new();
        for tree in &self.trees {
            for (var, imp) in tree.get_importance() {
                if tree_imps.contains_key(&var) {
                    let n_imp = &mut imp.to_vec();
                    let mut n_vec = tree_imps[&var].to_vec();
                    n_vec.append(n_imp);
                    tree_imps.remove(&var);
                    tree_imps.insert(var, n_vec);
                }
                else {
                    tree_imps.insert(var, imp);
                }
            }
        }
        for (var, imp) in tree_imps {
            let mean: f32 = imp.iter().sum::<f32>() / imp.len() as f32;
            println!("{:?}\t{:?}", var, mean);
        }
    }
}

fn make_tree(gm: &matrix::GenoMatrix, n_vars: &f64, n_subj: &f64, min_node_size: &i32, random_seed: &u64) -> tree::Node {
    let sample = gm.make_slice(n_vars, n_subj, random_seed);
    let data = gm.get_slice_data(&sample, random_seed);
    tree::Node::grow(data, *min_node_size, sample)
}

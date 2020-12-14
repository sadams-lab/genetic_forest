// Builds and analyzes random forest


use crate::tree;
use crate::matrix;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use std::collections::HashMap;


// A forest is simply a collection of trees (duh)
pub struct Forest {
    pub trees: Vec<tree::Node>
}

impl Forest {
    pub fn grow(gm: &matrix::GenoMatrix, n_tree: i32, mtry: f64, min_node_size: i32, subj_fraction: f64, threads: usize) -> Self {
        let var_sample_size: f64 = mtry * gm.n_genotypes;
        let subj_sample_size: f64 = subj_fraction * gm.n_subjects;
        let tp = ThreadPoolBuilder::new().num_threads(threads);
        match tp.build_global() {
            Ok(_) => eprintln!("Threads initialized successfully"),
            Err(e) => panic!("Error in threads initialization: {}", e),
        }
        let t: Vec<tree::Node> = (0..n_tree).into_par_iter()// Multithreaded piece, the tree building is recursive and can be independent
        .map(|_| -> tree::Node {
            make_tree(&gm, &var_sample_size, &subj_sample_size, &min_node_size)
        }).collect();
        Forest {trees: t}
    }

    pub fn re_grow(&mut self, gm: &matrix::GenoMatrix, n_tree: i32, mtry: f64, min_node_size: i32, subj_fraction: f64) {
        let var_sample_size: f64 = mtry * gm.n_genotypes;
        let subj_sample_size: f64 = subj_fraction * gm.n_subjects;
        let t: Vec<tree::Node> = (0..n_tree).into_par_iter()// Multithreaded piece, the tree building is recursive and can be independent
        .map(|_| -> tree::Node {
            make_tree(&gm, &var_sample_size, &subj_sample_size, &min_node_size)
        }).collect();
        self.trees = t;
    }

    fn importance(&self) -> HashMap<usize, Vec<f32>> {
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
        };
        tree_imps
    }

    pub fn mask_vars(&self, cutoff: f32) -> Vec<usize> {
        let mut vars: Vec<usize> = Vec::new();
        let tree_imps = self.importance();
        for (var, imp) in tree_imps {
            let mean: f32 = imp.iter().sum::<f32>() / imp.len() as f32;
            if mean < cutoff {
                vars.push(var);
            }
        }
        vars
    }

    pub fn pick_vars(&self, cutoff: f32) -> Vec<usize> {
        let mut vars: Vec<usize> = Vec::new();
        let tree_imps = self.importance();
        for (var, imp) in tree_imps {
            let mean: f32 = imp.iter().sum::<f32>() / imp.len() as f32;
            if mean > cutoff {
                vars.push(var);
            }
        }
        vars
    }
    
    pub fn print_var_importance(&self) {
        let tree_imps = self.importance();
        for (var, imp) in tree_imps {
            let mean: f32 = imp.iter().sum::<f32>() / imp.len() as f32;
            println!("{:?}\t{:?}", var, mean);
        }
    }
}

fn make_tree(gm: &matrix::GenoMatrix, n_vars: &f64, n_subj: &f64, min_node_size: &i32) -> tree::Node {
    // Connection to the tree lib for making the decision trees
    let sample = gm.make_slice(n_vars, n_subj);
    let data = gm.get_slice_data(&sample);
    tree::Node::grow(data, *min_node_size, sample)
}

// Builds and analyzes random forest

use crate::tree;
use crate::matrix;

use rayon::prelude::*;

use std::io;
use std::collections::HashMap;

pub struct HyperParameters {
    n_tree: i32,
    mtry: f64,
    min_node_size: i32,
    subj_fraction: f64
}
pub struct Forest {
    hyperparameters: HyperParameters,
    pub trees: Option<Vec<tree::Node>>
}

impl HyperParameters {

    pub fn set(n_tree: i32, mtry: f64, min_node_size: i32, subj_fraction: f64) -> Self {
        HyperParameters {
            n_tree: n_tree,
            mtry: mtry,
            min_node_size: min_node_size,
            subj_fraction: subj_fraction
        }
    }
}

impl Forest {

    pub fn new(hp: HyperParameters) -> Self {

        Forest {
            hyperparameters: hp,
            trees: None
        }
    }

    pub fn grow(&mut self, gm: &matrix::GenoMatrix) -> Result<(), io::Error> {
        let t: Vec<tree::Node> = (0..self.hyperparameters.n_tree).into_par_iter()
        .map(|_| -> tree::Node {
            match make_tree(gm, &self.hyperparameters) {
                Ok(tree) => return tree,
                Err(_) => return tree::Node::empty_node(),
            };
        }).collect();
        self.trees = Some(t);
        Ok(())
    }

    fn importance(&self) -> HashMap<usize, Vec<f32>> {
        let mut tree_imps: HashMap<usize, Vec<f32>> = HashMap::new();
        for tree in self.trees.as_ref().unwrap() {
            if !tree.is_empty {
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

fn make_tree(gm: &matrix::GenoMatrix, hp: &HyperParameters) -> Result<tree::Node, io::Error> {
    // Connection to the tree lib for making the decision trees
    let sample = gm.make_slice(hp.mtry, hp.subj_fraction);
    let data = gm.get_slice_data(&sample);
    Ok(tree::Node::grow(data, hp.min_node_size, sample))
}

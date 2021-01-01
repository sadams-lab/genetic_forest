// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

//! Random Forest Algorithm
//! Manages the creation and organization of decision trees

use crate::tree;
use crate::matrix;
use crate::variants;

use rayon::prelude::*;

use std::io;
use std::collections::HashMap;

/// Hyperparameters
/// n_tree = number of trees
/// mtry = fraction of variants to be selected for each tree
/// min_node_size = minimum number of samples in a node to be considered for a split
/// subj_fraction = fraction of subjects to be selected for each tree
pub struct HyperParameters {
    pub n_tree: i32,
    pub mtry: f64,
    pub max_depth: i32,
    pub subj_fraction: f64
}

pub struct Forest {
    hyperparameters: HyperParameters,
    pub trees: Option<Vec<tree::Node>>
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
    
    /// print the variant + importance to stdout
    pub fn print_var_importance(&self, variants: &Vec<variants::Variant>) {
        let tree_imps = self.importance();
        for (var, imp) in tree_imps {
            let mean: f32 = imp.iter().sum::<f32>() / imp.len() as f32;
            println!("{:?}\t{:?}", variants[var].id, mean);
        }
    }
}



/// Connection to the tree lib for making the decision trees
/// Outside of impl block since it is 'kind of' an independent operator
/// that spawns / returns the tree
fn make_tree(gm: &matrix::GenoMatrix, hp: &HyperParameters) -> Result<tree::Node, io::Error> {
    let sample = gm.make_slice(hp.mtry, hp.subj_fraction);
    let data = gm.get_slice_data(&sample);
    let tree_data = tree::NodeData {
        phenos: data.0,
        phenos_shuffle: data.1,
        genos: data.2
    };
    Ok(tree::Node::grow(tree_data, hp.max_depth, sample))
}

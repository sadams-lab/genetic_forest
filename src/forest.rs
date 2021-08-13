// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

//! Random Forest Algorithm
//! Manages the creation and organization of decision trees

use crate::tree;
use crate::matrix;
use crate::variants;
use crate::statistics;

use rayon::prelude::*;
use indicatif::ParallelProgressIterator;

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
    pub subj_fraction: f64,
    pub continuous_outcome: bool
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

    pub fn update_hyperparameters(&mut self, hp: HyperParameters) {
        self.hyperparameters = hp
    }

    pub fn grow(&mut self, gm: &matrix::GenoMatrix) -> Result<(), io::Error> {
        let t: Vec<tree::Node> = (0..self.hyperparameters.n_tree).into_par_iter()
        .progress_count(self.hyperparameters.n_tree as u64)
        .map(|_| -> tree::Node {
            match make_tree(gm, &self.hyperparameters) {
                Ok(tree) => return tree,
                Err(_) => return tree::Node::empty_node()
            };
        }).collect();
        self.trees = Some(t);
        Ok(())
    }

    fn get_var_importances(&self) -> HashMap<usize, f64> {
        let mut tree_imps: HashMap<usize, Vec<f64>> = HashMap::new();
        let mut var_imps: HashMap<usize, f64> = HashMap::new();
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
        for (var, imp) in tree_imps {
            var_imps.insert(var, statistics::mean(&imp.iter().map(|x| x).collect()));
        };
        var_imps
    }


    pub fn keep_vars(&self, p_keep: f64) -> Vec<usize> {
        let mut vars: Vec<usize> = Vec::new();
        let tree_imps = self.get_var_importances();
        let importances: Vec<&f64> = tree_imps.values().collect();
        let imp_mean: f64 = statistics::mean(&importances);
        let imp_sd: f64 = statistics::std_deviation(&importances);
        let cutoff: f64 = statistics::get_cutoff(&imp_sd, &imp_mean, &p_keep);
        for (var, imp) in tree_imps {
            if imp >= cutoff {
                vars.push(var);
            }
        }
        vars
    }
    
    /// print the variant + importance to stdout
    pub fn print_var_importance(&self, variants: &Vec<variants::Variant>) {
        let tree_imps = self.get_var_importances();
        for (var, imp) in tree_imps {
            println!("{}\t{:?}", variants[var].id, imp);
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
    Ok(tree::Node::grow(tree_data, hp.max_depth, sample, hp.continuous_outcome))
}

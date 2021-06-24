// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

//! Decision tree algorithm
//! Defines an abstract tree represented by nested nodes

use crate::matrix;
use crate::utils;
use crate::statistics;

use std::collections::HashMap;


/// Defines the node, which is the abstraction of the decision tree
/// Where the top level node is allocated to the stack, and subsequent
/// nodes on the heap. In order to 'draw' a tree, one must walk down recursive nodes
pub struct Node {
    pub score: f32, // gini if binary outcome, sd reduction if continuous
    pub is_empty: bool,
    pub n: usize, // Number of subjects in tree
    pub neg: bool,
    pub var: usize,
    pub node_n: usize, // Number of subjects in node, needed for importance calculation of parent node
    pub left: Option<Box<Node>>,
    pub right: Option<Box<Node>>
}

/// tree data
/// Contains data that is passed to create
/// a decision tree
pub struct NodeData<'a> {
    pub phenos: Vec<&'a f64>, // phenotype vector
    pub phenos_shuffle: Vec<&'a f64>, // vector of shuffled phenotypes (randomly, for impurity correction)
    pub genos: Vec<Vec<&'a u8>> // nested vector of genotypes
}

impl Node {

    /// The function that creates a tree
    pub fn grow(node_data: NodeData, max_depth: i32, ms: matrix::GenoMatrixSlice, continuous_outcome: bool) -> Self {
        let depth: i32 = 0;
        let node_n = node_data.phenos.len();
        return Node::new_node(&node_data, &ms, &max_depth, &depth, node_n, 1., &continuous_outcome);
    }

    /// An empty node, called internally to allow for terminal nodes that stop growth
    pub fn empty_node() -> Self {
        Node {
            score: std::f32::NAN,
            is_empty: true,
            n: 0,
            neg: true,
            var: 0,
            node_n: 0,
            left: None,
            right: None
        }
    }

    /// Print a tree to stdout by iterating over recursive nodes
    pub fn print(&self, above: &usize, side: &str) {
        println!("{:?}\t{:?}\t{}\t{:?}", self.var, above, side, self.score);
        match &self.left {
            Some(n) => n.print(&self.var, &"left"),
            None => ()
        }
        match &self.right {
            Some(n) => n.print(&self.var, &"right"),
            None => ()
        }
    }

    /// Calculate the importance of each variant used in a tree
    pub fn get_importance(&self) -> HashMap<usize, Vec<f32>> {
        // variant importances are stored in a hashmap with index: <importance, importance>
        // we store the importances as a vector because one feature might be selected multiple times in 
        // a single tree.
        let mut var_imp: HashMap<usize, Vec<f32>> = HashMap::new();
        /// Nested function to do the calculation
        fn imp(n: &Node, vi: &mut HashMap<usize, Vec<f32>>) {
            let mut importance: f32 = (n.node_n as f32 / n.n as f32) * (n.score);
            match &n.left {
                Some(nl) => {
                    if &nl.score > &0. {
                        importance -= (&nl.node_n / n.n) as f32 * (&nl.score);
                        imp(&nl, vi);
                    }
                },
                None => ()
            };
            match &n.right {
                Some(nr) => {
                    if &nr.score > &0. {
                        importance -= (&nr.node_n / n.n) as f32 * (&nr.score);
                        imp(&nr, vi);
                    }
                },
                None => ()
            };
            if n.neg {
                // was this calculated based on shuffled subjects?
                // If so, it is a negative contribution
                importance = importance * -1.;
            }
            if vi.contains_key(&n.var) {
                let mut n_varvec = vi[&n.var].to_vec();
                n_varvec.push(importance);
                vi.remove(&n.var);
                vi.insert(n.var, n_varvec);
            }
            else {
                vi.insert(n.var, vec![importance]);
            }
        }
        imp(self, &mut var_imp);
        return var_imp
    }

    /// Recursive function for building the tree
    /// Kicked off when a new tree is created
    fn new_node(node_data: &NodeData, ms: &matrix::GenoMatrixSlice, max_depth: &i32, depth: &i32, n: usize, parent_score: f32, continuous_outcome: &bool) -> Node {
        let new_depth = depth + 1;
        let mut best_score_index: usize = 0;
        let mut scores: Vec<f32> = Vec::new(); // Vector of per-genotype scores
        for k in &node_data.genos {
            // Each genotype vector has score calculated based on the actual phenos (score1)
            // and the shuffled phenotypes (score2)
            if *continuous_outcome {
                let score = calc_sdr(&node_data.phenos, k);
                let score2 = calc_sdr(&node_data.phenos_shuffle, k);
                match score {
                    // if shuffled pheno score is better than the actual pheno, then use that one
                    // but make it negative to indicate that it is to be a penalty rather than a contributor
                    // to overall importance
                    _x if score == 0. => (),
                    _x if score > score2 => {scores.push(score)},
                    _x if score <= score2 => {scores.push(score2 * -1.)},
                    _ => ()
                };
                best_score_index = utils::get_max_index(&scores);

            } else {
                let score = calc_gini(&node_data.phenos, k, &parent_score);
                let score2 = calc_gini(&node_data.phenos_shuffle, k, &parent_score);
                if score < score2 {
                    scores.push(score);
                }
                else { 
                    // if shuffled pheno score is better than the actual pheno, then use that one
                    // but make it negative to indicate that it is to be a penalty rather than a contributor
                    // to overall importance
                    scores.push(-1. * score2);
                }

                best_score_index = utils::get_min_index(&scores); // get lowest Gini score
            };
        }
        if scores.len() == 0 {
            return Node::empty_node(); 
        }
        let mut score = scores[best_score_index];
        let mut neg: bool = false;
        if score < 0. {
            neg = true;
            score = score * -1.;
        }
        let left_indices: &Vec<bool> = &node_data.genos[best_score_index].iter().map(|r| **r==0).collect();
        let right_indices: &Vec<bool> = &node_data.genos[best_score_index].iter().map(|r| **r==1).collect();
        let new_node_data = &node_data.split(&left_indices, &right_indices);
        if depth > max_depth {
            return Node {
                score: score,
                is_empty: false,
                n: n,
                neg: neg,
                node_n: node_data.phenos.len(),
                var: ms.genotype_ids[best_score_index],
                left: None,
                right: None
            };
        }
        Node {
            score: score,
            is_empty: false,
            n: n,
            neg: neg,
            node_n: node_data.phenos.len(),
            var: ms.genotype_ids[best_score_index],
            left: Some(Box::new(Node::new_node(&new_node_data.0, ms, max_depth, &new_depth, n, score, continuous_outcome))),
            right: Some(Box::new(Node::new_node(&new_node_data.1, ms, max_depth, &new_depth, n, score, continuous_outcome)))
        }
    }
}


/// Implementation of node data
/// Methods to manage the data that is used in the growing tree
impl<'a> NodeData<'a> {

    fn split(&self, left_indices: &Vec<bool>, right_indices: &Vec<bool>) -> (Self, Self) {
        let new_genos = self._split_genos(left_indices, right_indices);
        let new_phenos = self._split_phenos(left_indices, right_indices);
        let new_phenos_shuffle = self._split_phenos_shuffle(left_indices, right_indices);
        (
            NodeData {
                genos: new_genos.0,
                phenos: new_phenos.0,
                phenos_shuffle: new_phenos_shuffle.0

            }, 
            NodeData {
                genos: new_genos.1,
                phenos: new_phenos.1,
                phenos_shuffle: new_phenos_shuffle.1
            }
        )

    }

    fn _split_genos(&self, left_indices: &Vec<bool>, right_indices: &Vec<bool>) -> (Vec<Vec<&'a u8>>, Vec<Vec<&'a u8>>) {
        let mut left_genotypes: Vec<Vec<&u8>> = Vec::new();
        let mut right_genotypes: Vec<Vec<&u8>> = Vec::new();
        for g in &self.genos {
            let mut ind: usize = 0;
            let mut left_geno: Vec<&u8> = Vec::new();
            let mut right_geno: Vec<&u8> = Vec::new();
            for i in left_indices.iter().zip(right_indices.iter()) {
                if *i.0 {
                    left_geno.push(g[ind]);
                }
                if *i.1 {
                    right_geno.push(g[ind]);
                }
                ind += 1;
            }
            left_genotypes.push(left_geno);
            right_genotypes.push(right_geno)
        };
        return (left_genotypes, right_genotypes);
    }
    
    fn _split_phenos(&self, left_indices: &Vec<bool>, right_indices: &Vec<bool>) -> (Vec<&'a f64>, Vec<&'a f64>)  {
        let mut left_phenotypes: Vec<&f64> = Vec::new();
        let mut right_phenotypes: Vec<&f64> = Vec::new();
        let mut ind: usize = 0;
        for i in left_indices.iter().zip(right_indices.iter()) {
            if *i.0 {
                left_phenotypes.push(self.phenos[ind]);
            };
            if *i.1 {
                right_phenotypes.push(self.phenos[ind]);
            }
            ind += 1;
        };
        return (left_phenotypes, right_phenotypes);
    }

    fn _split_phenos_shuffle(&self, left_indices: &Vec<bool>, right_indices: &Vec<bool>) -> (Vec<&'a f64>, Vec<&'a f64>)  {
        let mut left_phenotypes: Vec<&f64> = Vec::new();
        let mut right_phenotypes: Vec<&f64> = Vec::new();
        let mut ind: usize = 0;
        for i in left_indices.iter().zip(right_indices.iter()) {
            if *i.0 {
                left_phenotypes.push(self.phenos_shuffle[ind]);
            };
            if *i.1 {
                right_phenotypes.push(self.phenos_shuffle[ind]);
            }
            ind += 1;
        };
        return (left_phenotypes, right_phenotypes);
    }
}

pub fn calc_sdr(p: &Vec<&f64>, g: &Vec<&u8>) -> f32 {
    let top_sd = statistics::std_deviation(p)as f32;
    let mut g0vec: Vec<&f64> = Vec::new();
    let mut g1vec: Vec<&f64> = Vec::new();
    for pg in p.iter().zip(g.iter()) {
        let (pi, gi) = pg;
        match *gi {
            0 => g0vec.push(pi),
            1 => g1vec.push(pi),
            _ => panic!("variable mismatch?, {}, {}", pi, gi),
        }
    }
    let g0sd = statistics::std_deviation(&g0vec) as f32;
    let g1sd = statistics::std_deviation(&g1vec) as f32;
    let sd_weighted = match g.len() {
        0 => 0.,
        _ => {(g0sd * (g0vec.len() / g.len()) as f32) + (g1sd * (g1vec.len() / g.len()) as f32)}
    };
    if sd_weighted > top_sd {
        return 0.
    }
    top_sd - sd_weighted

}

pub fn calc_gini(p: &Vec<&f64>, g: &Vec<&u8>, parent_score: &f32) -> f32 {
    /* 
    vector of genotypes (0,1,2) (g)
    vector of phenotypes (0,1) (p)

    geno 0 = pheno 0
    geno 0 = pheno 1
    ---
    get pheno where p = 0
    get pheno where p = 1

    */

    let pi: f32 = 2.;

    let mut p0g0: f32 = 0.;
    let mut p1g0: f32 = 0.;
    let mut p0g1: f32 = 0.;
    let mut p1g1: f32 = 0.;

    for pg in p.iter().zip(g.iter()) {
        let (pi, gi) = pg;
        match (**pi as u64, *gi) {
            (0, 0) => p0g0 += 1.,
            (1, 0) => p1g0 += 1.,
            (0, 1) => p0g1 += 1.,
            (1, 1) => p1g1 += 1.,
            _ => panic!("variable mismatch?, {}, {}", pi, gi),
        }
    }

    let p_inst: f32 = p.len() as f32;

    let g0_g = parent_score - (p0g0 / (p0g0+p1g0)).powf(pi) - (p1g0 / (p0g0+p1g0)).powf(pi);
    let g1_g = parent_score - (p0g1 / (p0g1+p1g1)).powf(pi) - (p1g1 / (p0g1+p1g1)).powf(pi);

    let gi = (((p0g0+p1g0)/p_inst) * g0_g) + (((p0g1+p1g1)/p_inst) * g1_g);
    if f32::is_nan(gi) {
        return 1.
    };
    return gi
}

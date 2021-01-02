//! Genetic algorithm for the optimization of a standard forest

use crate::forest;

use std::collections::HashMap;

/// Paired selection frequency
fn paired_selection_frequency(forest: &forest::Forest) {
    // Assemble a co-occurrence vector for each genotype
    let mut co_occurrence: Vec<Vec<i64>> = Vec::new();
    let mut variant_index: HashMap<usize, usize> = HashMap::new(); // so that it can map variant index with array index
    let variants = forest.get_vars();
    let mut i: usize = 0;
    for &v in variants {
        co_occurrence.push(Vec::new());
        variant_index.insert(v, i);
        i += 1;
    }
    for tree in forest.trees {
        // todo: best way to efficiently walk trees and populate this matrix
    }

}

/* paired_selection_freq <- function(var1, var2, ranger_obj) {
    var_names <- ranger_obj$forest$independent.variable.names
    n_trees <- final_model$forest$num.trees
    trees <- lapply(ranger_obj$forest$split.varIDs, FUN = function(x) {var_names[x]})
    var1_count <- sum(mapply(trees, FUN = function(x, y) {y %in% x}, y = var1 ))
    var2_count <- sum(mapply(trees, FUN = function(x, y) {y %in% x}, y = var2 ))
    var12_actual_on <- sum(mapply(trees, FUN = function(x, y, z) {y %in% x & z %in% x}, z = var1, y = var2 ))
    var12_prob_on <- (var1_count / n_trees) * (var2_count / n_trees)
    print((var12_actual_on/n_trees) - var12_prob_on)
    bi <- binom.test(c(var12_actual_on, n_trees - var12_actual_on), var12_prob_on, alternative = "greater")
    p <- bi$p.value
    return(p)
} */
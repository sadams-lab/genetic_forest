// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

pub mod reader;
pub mod matrix;
pub mod tree;
pub mod forest;
pub mod utils;
pub mod variants;
pub mod statistics;

use argparse::{ArgumentParser, StoreTrue, Store};


fn main() {
    let mut verbose: bool = false;
    let mut file_path: String = "".to_string();
    let mut variant_file_path: String = "".to_string();
    let mut n_tree: i32 = 10;
    let mut mtry: f64 = 0.333;
    let mut max_depth: i32 = 5;
    let mut subj_fraction: f64 = 0.666;
    let mut n_iter: usize = 1; 
    let mut var_cutoff: f32 = 0.;
    let mut continuous_outcome: bool = false;
    let mut output_forest: bool = false;
    let mut threads: usize = 1;
    {
        let mut ap = ArgumentParser::new();
        ap.set_description("Find genetic interactions from GWAS-scale data");
        ap.refer(&mut verbose)
        .add_option(&["-v", "--verbose"], StoreTrue,"Provide more verbose output");
        ap.refer(&mut file_path)
        .add_option(&["-f", "--file-path"], Store, "Path to input file.");
        ap.refer(&mut variant_file_path)
        .add_option(&["-x", "--variant-file-path"], Store, "Path to file with variant IDs");
        ap.refer(&mut n_tree)
        .add_option(&["-u", "--num-tree"], Store, "Number of trees to train.");
        ap.refer(&mut mtry)
        .add_option(&["-y", "--mtry-fraction"], Store, "MTRY fraction.");
        ap.refer(&mut max_depth)
        .add_option(&["-o", "--max-depth"], Store, "Max tree depth.");
        ap.refer(&mut subj_fraction)
        .add_option(&["-s", "--subject-fraction"], Store, "Fraction of subjects in each random sample.");
        ap.refer(&mut n_iter)
        .add_option(&["-r", "--iterations"], Store, "Number of iterations of the genetic algorithm.");
        ap.refer(&mut var_cutoff)
        .add_option(&["-c", "--cutoff"], Store, "Variant importance cutoff");
        ap.refer(&mut threads)
        .add_option(&["-t", "--threads"], Store, "Number of threads to use.");
        ap.refer(&mut output_forest)
        .add_option(&["--output-forest"], StoreTrue, "Output forest");
        ap.refer(&mut continuous_outcome)
        .add_option(&["--continuous-outcome"], StoreTrue, "Outcome is continuous (vs binary)");
        ap.parse_args_or_exit();
    }
    let filetype: u8 = input_file_type(&file_path);
    let mut data = match filetype {
        1 => reader::read_matrix_csv(&file_path, &",", &continuous_outcome),
        2 => reader::read_matrix_csv(&file_path, &"\t", &continuous_outcome),
        _ => panic!("Filetype not supported!"),
    };
    let variants = match filetype {
        1 => reader::read_variant_table(&variant_file_path, &",", &data.n_genotypes),
        2 => reader::read_variant_table(&variant_file_path, &"\t", &data.n_genotypes),
        _ => panic!("Filetype not supported!"),
    };
    utils::make_thread_pool(threads);
    eprintln!("Growing initial forest");
    let hp = forest::HyperParameters {
        n_tree: n_tree, 
        mtry: mtry, 
        max_depth: max_depth, 
        subj_fraction: subj_fraction,
        continuous_outcome: continuous_outcome
    };
    let mut f = forest::Forest::new(hp);
    match f.grow(&data) {
        Ok(_) => (),
        Err(err) => {
            println!("Error in initial forest growth: {}. Quitting now!", err);
            std::process::exit(1);
        }
    }
    let k_vars = f.keep_vars(var_cutoff);
    eprintln!("Keeping {:?} variants and initiating iterative grow and prune.", &k_vars.len());
    &data.set_genotype_indices(k_vars);
    for n in 1..n_iter {
        match f.grow(&data) {
            Ok(_) => (),
            Err(err) => {
                println!("Error in iteration {:?}: {}; breaking and returning results!", n, err);
                break
            }
        }
        let f_vars = f.mask_vars(var_cutoff);
        eprintln!("Masking {:?} variants and growing forest {:?} of {:?}.", &f_vars.len(), n, n_iter);
        &data.mask_vars(f_vars);
    }
    println!("#OVERALL_IMPORTANCE");
    f.print_var_importance(&variants);
    if output_forest {
        //Output trees
        for tree in f.trees.unwrap() {
            println!("#TREE");
            tree.print(&0, "0");
        };
    }
}

/// Determine input filetype based on suffix
fn input_file_type(filename: &str) -> u8 {
    let file_split: Vec<_> = filename.split(".").map(|s| s).collect();
    let suffix = file_split[file_split.len() - 1];
    match suffix {
        "csv" => 1,
        "tsv" => 2,
        "gz" => 9,
        _ => 0,
    }
}

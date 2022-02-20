// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

pub mod reader;
pub mod matrix;
pub mod tree;
pub mod forest;
pub mod utils;
pub mod variants;
pub mod statistics;

use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about = "Find genetic interactions from GWAS-scale data", long_about = None)]
struct Args {
    #[clap(long, help="Provide more verbose output")]
    verbose: bool,
    #[clap(short, long, help="Path to input file")]
    file_path: String,
    #[clap(long, help="Path to file with variant names.")]
    variant_file_path: String,
    #[clap(long, help="Number of trees in selection forest.")]
    n_tree: i32,
    #[clap(long, help="MTRY fraction for selection forest.")]
    mtry: f64,
    #[clap(long, help="Max depth for selection forest.")]
    max_depth: i32,
    #[clap(long, help="Subject fraction for selection forest.")]
    subj_fraction: f64,
    #[clap(long, help="Number of trees for iterative forest.")]
    n_tree_2: i32,
    #[clap(long, help="MTRY fraction for iterative forest.")]
    mtry_2: f64,
    #[clap(long, help="Max depth for iterative forest.")]
    max_depth_2: i32,
    #[clap(long, help="Subject fraction for iterative forest.")]
    subj_fraction_2: f64,
    #[clap(long, help="Number of iterations for iterative forest.")]
    n_iter: usize,
    #[clap(long, help="Z score to keep variants after selection forest.")]
    z_keep: f64,
    #[clap(long, help="Outcome is a continuous variable.")]
    continuous_outcome: bool,
    #[clap(long, help="Write forest to stdout (very verbose output)")]
    output_forest: bool,
    #[clap(short, long, help="Number of threads to use.")]
    threads: usize
}


fn main() {
    let args = Args::parse();
    let filetype: u8 = input_file_type(&args.file_path);
    let mut data = match filetype {
        1 => reader::read_matrix_csv(&args.file_path, &",", &args.continuous_outcome),
        2 => reader::read_matrix_csv(&args.file_path, &"\t", &args.continuous_outcome),
        _ => panic!("Filetype not supported!"),
    };
    let variants = match filetype {
        1 => reader::read_variant_table(&args.variant_file_path, &",", &data.n_genotypes),
        2 => reader::read_variant_table(&args.variant_file_path, &"\t", &data.n_genotypes),
        _ => panic!("Filetype not supported!"),
    };
    utils::make_thread_pool(args.threads);
    eprintln!("Growing initial forest");
    let hp = forest::HyperParameters {
        n_tree: args.n_tree, 
        mtry: args.mtry, 
        max_depth: args.max_depth, 
        subj_fraction: args.subj_fraction,
        continuous_outcome: args.continuous_outcome
    };
    let hp2 = forest::HyperParameters {
        n_tree: args.n_tree_2, 
        mtry: args.mtry_2, 
        max_depth: args.max_depth_2, 
        subj_fraction: args.subj_fraction_2,
        continuous_outcome: args.continuous_outcome
    };
    let mut f = forest::Forest::new(hp);
    match f.grow(&data) {
        Ok(_) => (),
        Err(err) => {
            println!("Error in initial forest growth: {}. Quitting now!", err);
            std::process::exit(1);
        }
    }
    f.print_var_importance(&variants);
    let k_vars = f.keep_vars(args.z_keep);
    eprintln!("Keeping {:?} variants and initiating iterative grow and prune.", &k_vars.len());
    &data.set_genotype_indices(k_vars);
    f.update_hyperparameters(hp2);
    for n in 1..&args.n_iter + 1 {
        eprintln!("Growing forest {:?} of {:?}", n, &args.n_iter);
        match f.grow(&data) {
            Ok(_) => (),
            Err(err) => {
                println!("Error in iteration {:?}: {}; breaking and returning results!", n, err);
                break
            }
        }
        println!("## ITERATION: {}", n);
        println!("#OVERALL_IMPORTANCE");
        f.print_var_importance(&variants);
        if args.output_forest {
            //Output trees
            for tree in f.trees.as_ref().unwrap() {
                println!("#TREE");
                tree.print(&0, "0");
            };
        }
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

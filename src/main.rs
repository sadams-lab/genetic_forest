pub mod reader;
pub mod matrix;
pub mod tree;
pub mod forest;

use argparse::{ArgumentParser, StoreTrue, Store};


fn main() {
    let mut verbose: bool = false;
    let mut file_path: String = "".to_string();
    let mut n_tree: i32 = 10;
    let mut mtry: f64 = 0.333;
    let mut min_node_size: i32 = 2;
    let mut subj_fraction: f64 = 0.666;
    let mut threads: usize = 1;
    {
        let mut ap = ArgumentParser::new();
        ap.set_description("Find genetic interactions from GWAS-scale data");
        ap.refer(&mut verbose)
        .add_option(&["-v", "--verbose"], StoreTrue,"Provide more verbose output");
        ap.refer(&mut file_path)
        .add_option(&["-f", "--file-path"], Store, "Path to input file.");
        ap.refer(&mut n_tree)
        .add_option(&["-u", "--num-tree"], Store, "Number of trees to train.");
        ap.refer(&mut mtry)
        .add_option(&["-y", "--mtry-fraction"], Store, "MTRY fraction.");
        ap.refer(&mut min_node_size)
        .add_option(&["-o", "--min-node-size"], Store, "Minimum number of subjects in a node (surrogate for max-depth).");
        ap.refer(&mut subj_fraction)
        .add_option(&["-s", "--subject-fraction"], Store, "Fraction of subjects in each random sample.");
        ap.refer(&mut threads)
        .add_option(&["-t", "--threads"], Store, "Number of threads to use.");
        ap.parse_args_or_exit();
    }
    let filetype: u8 = input_file_type(&file_path);
    let data = match filetype {
        1 => reader::read_csv(&file_path, &","),
        2 => reader::read_csv(&file_path, &"\t"),
        _ => panic!("Filetype not supported!"),
    };
    let f = forest::Forest::grow(data, n_tree, mtry, min_node_size, subj_fraction, threads);
    f.var_importance();
}

fn input_file_type(filename: &str) -> u8 {
    // Determine if filetype is csv, tsv, or gz
    // TODO: add support for gz files
    let file_split: Vec<_> = filename.split(".").map(|s| s).collect();
    let suffix = file_split[file_split.len() - 1];
    match suffix {
        "csv" => 1,
        "tsv" => 2,
        "gz" => 9,
        _ => 0,
    }
}


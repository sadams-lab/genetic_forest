pub mod reader;
pub mod matrix;
pub mod tree;
pub mod forest;

use std::env;


fn main() {
    let args: Vec<String> = env::args().collect();
    let filetype: u8 = input_file_type(&args[1]);
    let data = match filetype {
        1 => reader::read_csv(&args[1], &","),
        2 => reader::read_csv(&args[1], &"\t"),
        _ => panic!("Input filetype not implemented!"),
    };
    let f = forest::Forest::grow(data, 10000, 0.01, 5, 0.333, 4);
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


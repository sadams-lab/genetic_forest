mod reader;
mod forest;

use std::env;


fn main() {
    let args: Vec<String> = env::args().collect();
    let filetype: u8 = input_file_type(&args[1]);
    println!("{:?}, {:?}", args, filetype);
    let data = match filetype {
        1 => reader::read_csv(&args[1]),
        _ => panic!("Input filetype not implemented!"),
    };
    let sample = data.random_select(&500., &250.);
    let phenos = &sample.0;
    for k in &sample.1 {
        let g = forest::gini(phenos, k);
        println!("{}", g);
    }
}

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

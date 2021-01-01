// simplest use case
// expects a csv file
// ID,PHENO,COV1,COV2,COV3

use crate::matrix;
use crate::variants;
use crate::utils;

use std::str;
use std::fs::File;

pub fn read_matrix_csv(path: &str, sep: &str) -> matrix::GenoMatrix {
    let mut beginning: csv::Position = csv::Position::new();
    beginning.set_line(1);
    let mut reader = make_reader(path, sep);
    let shape = utils::get_mat_size(&mut reader);
    match reader.seek(beginning) {
        Ok(_) => (),
        Err(e) => panic!("Error in file seek {:?}", e)
    }
    return matrix::GenoMatrix::new(&mut reader, shape);
}

pub fn read_variant_table(path: &str, sep: &str) -> Vec<variants::Variant> {
    let mut variants: Vec<variants::Variant> = Vec::new();
    let mut reader = make_reader(path, sep);
    for result in reader.records() {
        let record = result.unwrap();
        variants.push(variants::Variant::new(String::from(&record[1])));
    }
    variants
}

fn make_reader(path: &str, sep: &str) -> csv::Reader<File> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(why) => panic!("Something happened with the file: {}", why),
    };
    return csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(sep.as_bytes()[0])
        .comment(Some(b'#'))
        .from_reader(file);
}
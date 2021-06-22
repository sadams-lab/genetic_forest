// Copyright 2020 Solomon M. Adams, PharmD, PhD
// Licensed under the MIT license

//! Functions related to file reading

use crate::matrix;
use crate::variants;
use crate::utils;

use std::str;
use std::fs::File;

pub fn read_matrix_csv(path: &str, sep: &str) -> matrix::GenoMatrix {
    let mut beginning: csv::Position = csv::Position::new();
    beginning.set_line(1);
    let mut reader = match make_reader(path, sep) {
        Ok(reader) => reader,
        Err(why) => panic!(why)
    };
    let shape = utils::get_mat_size(&mut reader);
    match reader.seek(beginning) {
        Ok(_) => (),
        Err(e) => panic!("Error in file seek {:?}", e)
    }
    return matrix::GenoMatrix::new(&mut reader, shape);
}

pub fn read_variant_table(path: &str, sep: &str, n_genotypes: &f64) -> Vec<variants::Variant> {
    let mut variants: Vec<variants::Variant> = Vec::new();
    match make_reader(path, sep) {
        Ok(mut reader) => build_variant_array(&mut variants, &mut reader),
        Err(_) => build_dummy_variant_array(&mut variants, n_genotypes)
    };
    variants
}

fn build_dummy_variant_array(variants: &mut Vec<variants::Variant>, n_genotypes: &f64) {
    // When a variant id file is not provided
    for i in 0..n_genotypes.round() as i64 {
        variants.push(variants::Variant::new(i.to_string()));
    }
}
    
fn build_variant_array(variants: &mut Vec<variants::Variant>, reader: &mut csv::Reader<File>) {
    for result in reader.records() {
        let record = result.unwrap();
        variants.push(variants::Variant::new(String::from(&record[1])));
    }
}

fn make_reader(path: &str, sep: &str) -> Result<csv::Reader<File>, &'static str> {
    match File::open(path) {
        Ok(f) => Ok(csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(sep.as_bytes()[0])
            .comment(Some(b'#'))
            .from_reader(f)),
        Err(_) => Err("Unable to load file"),
    }
}

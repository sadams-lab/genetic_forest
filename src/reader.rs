// simplest use case
// expects a csv file
// ID,PHENO,COV1,COV2,COV3

mod matrix;

use std::str;
use std::mem;

pub fn read_csv(path: &str) -> i32 {
    let foo = matrix::GenoMatrix::new(path);
    1
}

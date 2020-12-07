// simplest use case
// expects a csv file
// ID,PHENO,COV1,COV2,COV3

use crate::matrix;

use std::str;

pub fn read_csv(path: &str) -> matrix::GenoMatrix {
    return matrix::GenoMatrix::new(path);
}

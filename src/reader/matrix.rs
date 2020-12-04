use std::fs::File;
use std::str;
use sprs::{CsMat, Shape, TriMat};

pub struct GenoMatrix {
    pub ids: Vec<String>,
    pub phenotypes: Vec<u8>,
    genotypes: CsMat<u8>
}

impl GenoMatrix {
    
    pub fn new(path: &str) -> Self {
        let mut row_ids: Vec<String> = Vec::new();
        let mut phenotypes: Vec<u8> = Vec::new();
        let mat_size: Shape = get_mat_size(path);
        let mut geno_mat = TriMat::new(mat_size);
        let mut rdr = make_reader(path);
        let mut rownum = 0;
        for result in rdr.records() {
            let mut colnum = 2;
            let record = result.unwrap();
            row_ids.push(record[0].to_string());
            phenotypes.push(record[1].parse::<u8>().unwrap());
            loop {
                if colnum == mat_size.1 + 2 {
                    break
                }
                geno_mat.add_triplet(rownum, colnum - 2, record[colnum].parse::<u8>().unwrap());
                colnum += 1;
            }
            rownum += 1
        }
        GenoMatrix{
            ids: row_ids, 
            phenotypes: phenotypes,
            genotypes: geno_mat.to_csr(),
        }
    }
}

fn make_reader(path: &str) -> csv::Reader<File> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(why) => panic!("Something happened with the file: {}", why),
    };
    return csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b',')
        .comment(Some(b'#'))
        .from_reader(file);
}

fn get_mat_size(path: &str) -> Shape {
    let mut ncols: usize = 0;
    let mut nrows: usize = 0;
    let mut rdr = make_reader(path);
    for result in rdr.byte_records() {
            let record = result.unwrap();
            if ncols == 0 {
                ncols = &record.len() - 2;
            };
            nrows += 1;
    }
    (nrows, ncols)
}

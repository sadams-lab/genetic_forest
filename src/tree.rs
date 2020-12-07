// Make a decision tree

use crate::matrix;

struct Node {
    gini: f32,
    var: usize,
    left: Box<Node>,
    right: Box<Node>
}

impl Node {
    pub fn grow(data: (Vec<&u8>,  Vec<Vec<&u8>>), max_depth: i32, ms: matrix::GenoMatrixSlice) -> Self {
        return new_node(data, &ms);
    }
}

fn new_node(data: (Vec<&u8>,  Vec<Vec<&u8>>), ms: &matrix::GenoMatrixSlice) -> Node {
    let mut ginis: Vec<f32> = Vec::new();
    let phenos = &data.0;
    for k in &data.1 {
        ginis.push(calc_gini(phenos, k));
    };
    let min_gini_index: usize = min_gini(&ginis);
    let gini = ginis[min_gini_index];
    let left_indices = data.1[min_gini_index].iter().map(|r| r==&&0).collect();
    let right_indices = data.1[min_gini_index].iter().map(|r| r==&&0).collect();
    let left_genotypes = genotypes_split(&data.1, &left_indices);
    let right_genotypes = genotypes_split(&data.1, &right_indices);
    let left_phenotypes = phenotypes_split(&data.0, &left_indices);
    let right_phenotypes = phenotypes_split(&data.0, &right_indices);
    Node {
        gini: gini,
        var: ms.genotype_ids[min_gini_index],
        left: Box::new(new_node((left_phenotypes, left_genotypes), ms)),
        right: Box::new(new_node((right_phenotypes, right_genotypes), ms))
    }
}

fn genotypes_split<'a>(genotypes: &Vec<Vec<&'a u8>>, indices: &Vec<bool>) -> Vec<Vec<&'a u8>> {
    let mut new_genotypes: Vec<Vec<&u8>> = Vec::new();
    for g in genotypes {
        let mut ind: usize = 0;
        let mut new_geno: Vec<&u8> = Vec::new();
        for i in indices {
            if *i {
                new_geno.push(g[ind]);
            }
        }
        new_genotypes.push(new_geno);
        ind += 1;
    };
    return new_genotypes;
}

fn phenotypes_split<'a>(phenotypes: &Vec<&'a u8>, indices: &Vec<bool>) -> Vec<&'a u8> {
    let mut new_phenotypes: Vec<&u8> = Vec::new();
    let mut ind: usize = 0;
    for i in indices {
        if *i {
            new_phenotypes.push(phenotypes[ind]);
        }
    ind += 1;
    };
    return new_phenotypes;
}


fn min_gini(ginis: &Vec<f32>) -> usize {
    // Take vector of floats, return the index of the minumum
    let mut min_gin: &f32 = &1.;
    let mut min_i: usize = 0;
    let mut i: usize = 0;
    for g in ginis {
        if g < min_gin {
            min_gin = g;
            min_i = i;
        };
        i += 1;
    }
    return min_i
}

pub fn calc_gini(p: &Vec<&u8>, g: &Vec<&u8>) -> f32 {
    /* 
    vector of genotypes (0,1,2) (g)
    vector of phenotypes (0,1) (p)

    geno 0 = pheno 0
    geno 0 = pheno 1
    ---
    get pheno where p = 0
    get pheno where p = 1

    */

    let pi: f32 = 2.;

    let mut p0g0: f32 = 0.;
    let mut p1g0: f32 = 0.;
    let mut p0g1: f32 = 0.;
    let mut p1g1: f32 = 0.;
    let mut p0g2: f32 = 0.;
    let mut p1g2: f32 = 0.;

    for pg in p.iter().zip(g.iter()) {
        let (pi, gi) = pg;
        match (*pi, *gi) {
            (0, 0) => p0g0 += 1.,
            (1, 0) => p1g0 += 1.,
            (0, 1) => p0g1 += 1.,
            (1, 1) => p1g1 += 1.,
            _ => panic!("variable mismatch?, {}, {}", pi, gi),
        }
    }

    let p_inst: f32 = p.len() as f32;

    let g0_g = 1. - (p0g0 / (p0g0+p1g0)).powf(pi) - (p1g0 / (p0g0+p1g0)).powf(pi);
    let g1_g = 1. - (p0g1 / (p0g1+p1g1)).powf(pi) - (p1g1 / (p0g1+p1g1)).powf(pi);

    let gi = (((p0g0+p1g0)/p_inst) * g0_g) + (((p0g1+p1g1)/p_inst) * g1_g);
    if f32::is_nan(gi) {
        return 0.
    };
    return gi

}
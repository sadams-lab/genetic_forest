// Make a decision tree
//TODO: make the shuffle/calc happen during the gini calculation!
use crate::matrix;

use std::collections::HashMap;

#[derive(Default)]
pub struct Node {
    pub gini: f32,
    pub var: usize,
    pub left: Option<Box<Node>>,
    pub right: Option<Box<Node>>
}

impl Node {
    pub fn grow(data: (Vec<&u8>, Vec<&u8>, Vec<Vec<&u8>>), min_count: i32, ms: matrix::GenoMatrixSlice) -> Self {
        return new_node(data, &ms, &min_count);
    }

    pub fn print(&self, above: &usize, side: &str) {
        println!("{:?}\t{}\t{:?}\t{:?}", above, side, self.var, self.gini);
        match &self.left {
            Some(n) => n.print(&self.var, &"left"),
            None => ()
        }
        match &self.right {
            Some(n) => n.print(&self.var, &"right"),
            None => ()
        }
    }

    pub fn get_importance(&self) -> HashMap<usize, Vec<f32>> {
        let mut var_imp: HashMap<usize, Vec<f32>> = HashMap::new();
        fn imp(n: &Node, vi: &mut HashMap<usize, Vec<f32>>) {
            if vi.contains_key(&n.var) {
                let mut n_varvec = vi[&n.var].to_vec();
                n_varvec.push(n.gini);
                vi.remove(&n.var);
                vi.insert(n.var, n_varvec);
            }
            else {
                vi.insert(n.var, vec![n.gini]);
            }
            match &n.left {
                Some(nl) => imp(&nl, vi),
                None => ()
            };
            match &n.right {
                Some(nr) => imp(&nr, vi),
                None => ()
            };
        }
        imp(self, &mut var_imp);
        return var_imp
    }
}


fn new_node(data: (Vec<&u8>, Vec<&u8>, Vec<Vec<&u8>>), ms: &matrix::GenoMatrixSlice, min_count: &i32) -> Node {
    let mut ginis: Vec<f32> = Vec::new();
    let phenos = &data.0;
    let phenos2 = &data.1;
    for k in &data.2 {
        let gini = calc_gini(phenos, k);
        let gini2 = calc_gini(phenos2, k);
        if gini < gini2 {
            ginis.push(gini);
        }
        else {
            ginis.push(-1. * gini2);
        }
    };
    let min_gini_index: usize = min_gini(&ginis);
    let gini = ginis[min_gini_index];
    let left_indices: &Vec<bool> = &data.2[min_gini_index].iter().map(|r| **r==0).collect();
    let right_indices: &Vec<bool> = &data.2[min_gini_index].iter().map(|r| **r==1).collect();
    let new_genotypes = genotypes_split(&data.2, &left_indices, &right_indices);
    let new_phenotypes = phenotypes_split(&data.0, &left_indices, &right_indices);
    let new_phenotypes_2 = phenotypes_split(&data.1, &left_indices, &right_indices);
    if sum_bool_vec(&left_indices) <= *min_count && sum_bool_vec(&right_indices) <= *min_count {
        return Node {
            gini: gini,
            var: ms.genotype_ids[min_gini_index],
            left: None,
            right: None
        };
    }
    else if sum_bool_vec(&left_indices) <= *min_count {
        return Node {
            gini: gini,
            var: ms.genotype_ids[min_gini_index],
            left: None,
            right: Some(Box::new(new_node((new_phenotypes.1, new_phenotypes_2.1, new_genotypes.1), ms, min_count)))
        };

    }
    else if sum_bool_vec(&right_indices) <= *min_count {
        return Node {
            gini: gini,
            var: ms.genotype_ids[min_gini_index],
            left: Some(Box::new(new_node((new_phenotypes.0, new_phenotypes_2.0, new_genotypes.0), ms, min_count))),
            right: None
        };
        
    }
    Node {
        gini: gini,
        var: ms.genotype_ids[min_gini_index],
        left: Some(Box::new(new_node((new_phenotypes.0, new_phenotypes_2.0, new_genotypes.0), ms, min_count))),
        right: Some(Box::new(new_node((new_phenotypes.1, new_phenotypes_2.1, new_genotypes.1), ms, min_count)))
    }
}

fn sum_bool_vec(v: &Vec<bool>) -> i32 {
    let mut s: i32 = 0;
    for b in v{
        if *b {
            s += 1;
        };
    }
    return s;
}

fn genotypes_split<'a>(genotypes: &Vec<Vec<&'a u8>>, left_indices: &Vec<bool>, right_indices: &Vec<bool>) -> (Vec<Vec<&'a u8>>, Vec<Vec<&'a u8>>) {
    let mut left_genotypes: Vec<Vec<&u8>> = Vec::new();
    let mut right_genotypes: Vec<Vec<&u8>> = Vec::new();
    for g in genotypes {
        let mut ind: usize = 0;
        let mut left_geno: Vec<&u8> = Vec::new();
        let mut right_geno: Vec<&u8> = Vec::new();
        for i in left_indices.iter().zip(right_indices.iter()) {
            if *i.0 {
                left_geno.push(g[ind]);
            }
            if *i.1 {
                right_geno.push(g[ind]);
            }
            ind += 1;
        }
        left_genotypes.push(left_geno);
        right_genotypes.push(right_geno)
    };
    return (left_genotypes, right_genotypes);
}

fn phenotypes_split<'a>(phenotypes: &Vec<&'a u8>, left_indices: &Vec<bool>, right_indices: &Vec<bool>) -> (Vec<&'a u8>, Vec<&'a u8>)  {
    let mut left_phenotypes: Vec<&u8> = Vec::new();
    let mut right_phenotypes: Vec<&u8> = Vec::new();
    let mut ind: usize = 0;
    for i in left_indices.iter().zip(right_indices.iter()) {
        if *i.0 {
            left_phenotypes.push(phenotypes[ind]);
        };
        if *i.1 {
            right_phenotypes.push(phenotypes[ind]);
        }
        ind += 1;
    };
    return (left_phenotypes, right_phenotypes);
}


fn min_gini(ginis: &Vec<f32>) -> usize {
    // Take vector of floats, return the index of the minumum
    let mut min_gin: f32 = 1.;
    let mut min_i: usize = 0;
    let mut i: usize = 0;
    for g in ginis {
        let abs_gini = g.abs();
        if abs_gini < min_gin {
            min_gin = *g;
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
        return 1.
    };
    return gi
}
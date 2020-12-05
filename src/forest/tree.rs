// Make a decision tree

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
            (0, 2) => p0g2 += 1.,
            (1, 2) => p1g2 += 1.,
            _ => panic!("variable mismatch?, {}, {}", pi, gi),
        }
    }

    let p_inst: f32 = p.len() as f32;

    let g0_g = 1. - (p0g0 / (p0g0+p1g0)).powf(pi) - (p1g0 / (p0g0+p1g0)).powf(pi);
    let g1_g = 1. - (p0g1 / (p0g1+p1g1)).powf(pi) - (p1g1 / (p0g1+p1g1)).powf(pi);
    let g2_g = 1. - (p0g2 / (p0g2+p1g2)).powf(pi) - (p1g2 / (p0g2+p1g2)).powf(pi);

    let gi = (((p0g0+p1g0)/p_inst) * g0_g) + (((p0g1+p1g1)/p_inst) * g1_g) + (((p0g2+p1g2)/p_inst) * g2_g);
    if f32::is_nan(gi) {
        return 0.
    };
    return gi

}
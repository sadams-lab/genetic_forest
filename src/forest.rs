mod tree;

pub fn gini(p: &Vec<&u8>, g: &Vec<&u8>) -> f32 {
    return tree::calc_gini(p, g);
}
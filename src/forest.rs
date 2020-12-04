mod tree;

pub fn foo() {
    let g: Vec<u8> = vec![0,1,2,1,1,2,1,0,1,0,1,2,1,2,2,1,1,0];
    let p: Vec<u8> = vec![0,1,0,1,1,0,1,0,1,1,1,0,0,1,0,1,0,0];
    let gini: f32 = tree::calc_gini(&p, &g);
    println!("Gini is: {}", gini);
}
//! Miscellaneious utilities shared throughout the genetic forest module

use rayon::ThreadPoolBuilder;
use sprs::Shape;

use std::fs::File;

/// Initialize the global thread pool, only needs to be called once per program run
/// No return, just needs to be called with an argument corresponding to the number of threads
/// Panics when something goes wrong here with (hopefully) a helpful error message
pub fn make_thread_pool(n_threads: usize) {
    // if n_threads is bigger than available threads, then rayon handles the checking and will just assign the max number
    // available in the thread pool
    let tp = ThreadPoolBuilder::new().num_threads(n_threads);
    match tp.build_global() {
        Ok(_) => eprintln!("Threads initialized successfully"),
        Err(e) => panic!("Error in threads initialization: {}", e),
    }
}

/// Get the sum of a vector of booleans
/// where true == 1 and false == 0
pub fn sum_bool_vec(v: &Vec<bool>) -> i32 {
    let mut s: i32 = 0;
    for b in v{
        if *b {
            s += 1;
        };
    }
    return s;
}

/// Take vector of floats, return the index of the minumum
pub fn get_min_index(vals: &Vec<f32>) -> usize {
    let mut min_val: f32 = 0.;
    let mut min_i: usize = 0;
    let mut i: usize = 0;
    for g in vals {
        let abs_val = g.abs(); // get absolute value since they might be negative
        if abs_val < min_val {
            min_val = *g;
            min_i = i;
        };
        i += 1;
    }
    return min_i
}

/// Take vector of floats, return the index of the maximum
pub fn get_max_index(vals: &Vec<f32>) -> usize {
    let mut max_val: f32 = 0.;
    let mut max_i: usize = 0;
    let mut i: usize = 0;
    for g in vals {
        let abs_val = g.abs(); // get absolute value since they might be negative
        if abs_val > max_val {
            max_val = *g;
            max_i = i;
        };
        i += 1;
    }
    return max_i
}

/// Get the size of a reader object
/// In number of lines
pub fn get_mat_size(rdr: &mut csv::Reader<File>) -> Shape {
    let mut ncols: usize = 0;
    let mut nrows: usize = 0;
    for result in rdr.byte_records() {
            let record = result.unwrap();
            if ncols == 0 {
                ncols = &record.len() - 2;
            };
            nrows += 1;
    }
    (nrows, ncols)
}
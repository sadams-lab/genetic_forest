//! Miscellaneious utilities shared throughout the genetic forest module

use rayon::ThreadPoolBuilder;

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

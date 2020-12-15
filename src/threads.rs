use rayon::ThreadPoolBuilder;

pub fn make_thread_pool(n_threads: usize) {
    let tp = ThreadPoolBuilder::new().num_threads(n_threads);
    match tp.build_global() {
        Ok(_) => eprintln!("Threads initialized successfully"),
        Err(e) => panic!("Error in threads initialization: {}", e),
    }
}

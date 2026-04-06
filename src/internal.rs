use std::time::Duration;

use simple_sds::binaries;

//-----------------------------------------------------------------------------

pub fn report_results(queries: usize, total_len: usize, total_occs: usize, duration: Duration) {
    let us = (duration.as_micros() as f64) / (queries as f64);
    let ns = (duration.as_nanos() as f64) / (total_len as f64);
    let occs = (total_occs as f64) / (queries as f64);
    eprintln!("Time:        {:.3} seconds ({:.3} us/query, {:.1} ns/node)", duration.as_secs_f64(), us, ns);
    eprintln!("Occurrences: {} total ({:.3} per query)", total_occs, occs);
    eprintln!();
}

pub fn report_memory_usage() {
    match binaries::peak_memory_usage() {
        Ok(bytes) => {
            let (size, unit) = binaries::human_readable_size(bytes);
            eprintln!("Peak memory usage: {:.3} {}", size, unit);
        },
        Err(f) => {
            eprintln!("{}", f);
        },
    }
}

//-----------------------------------------------------------------------------

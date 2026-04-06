#![allow(
    clippy::uninlined_format_args,
    clippy::new_without_default
)]

use gbz::{GBWT, GBWTBuilder};
use gbz::FullPathName;
use gbz::internal;

use simple_sds::binaries;
use simple_sds::serialize::{self, Serialize};

use std::time::Instant;
use std::{env, process};

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start = Instant::now();
    let config = Config::new();

    let filename = config.filename.as_ref().unwrap();
    eprintln!("Loading GBWT index {}", filename);
    let truth: GBWT = serialize::load_from(filename).map_err(|e| format!("Failed to load GBWT index: {}", e))?;
    let (size, units) = binaries::human_readable_size(truth.size_in_bytes());
    eprintln!("{} sequences of total length {}; index size {:.3} {}", truth.sequences(), truth.len(), size, units);

    // Construction.
    let construction_start = Instant::now();
    let mut builder = GBWTBuilder::new(truth.is_bidirectional(), truth.has_metadata(), config.buffer_size);
    let increment = if truth.is_bidirectional() { 2 } else { 1 };
    let mut sequence_id = 0;
    while sequence_id < truth.sequences() {
        let sequence = truth.sequence(sequence_id).ok_or_else(||
            format!("Failed to extract sequence id {}", sequence_id)
        )?;
        let sequence: Vec<usize> = sequence.collect();
        let path_name = if let Some(metadata) = truth.metadata() {
            FullPathName::from_metadata(metadata, sequence_id / increment)
        } else {
            None
        };
        builder.insert(&sequence, path_name)?;
        sequence_id += increment;        
    }
    let index = builder.build()?;
    let seconds = construction_start.elapsed().as_secs_f64();
    let nodes_second = (truth.len() as f64) / seconds;
    eprintln!("Built GBWT index in {:.3} seconds ({:.3} nodes/second)", seconds, nodes_second);

    gbz::gbwt::builder::compare_gbwts(&index, &truth, "rebuild-gbwt")?;

    eprintln!("Used {:.3} seconds", start.elapsed().as_secs_f64());
    internal::report_memory_usage();
    eprintln!();
    Ok(())
}

//-----------------------------------------------------------------------------

pub struct Config {
    pub filename: Option<String>,
    pub buffer_size: usize,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("b", "buffer", "buffer size in nodes (default: 100M)", "INT");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        let mut config = Config {
            filename: None,
            buffer_size: 100_000_000
        };

        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] index.gbwt", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(arg) = matches.opt_str("b") {
            config.buffer_size = match binaries::parse_unsigned(&arg) {
                Ok(size) => size,
                Err(e) => {
                    eprintln!("--buffer: {}", e);
                    process::exit(1);
                }
            };
        }

        if let Some(filename) = matches.free.first() {
            config.filename = Some(filename.clone());
            config
        } else {
            let header = format!("Usage: {} [options] index.gbwt", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        }
    }
}

//-----------------------------------------------------------------------------

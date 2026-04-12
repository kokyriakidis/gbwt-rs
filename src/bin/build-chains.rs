#![allow(
    clippy::uninlined_format_args,
    clippy::new_without_default
)]

use gbz::GBZ;
use gbz::support::Chains;
use gbz::{algorithms, internal, support};

use simple_sds::serialize;

use std::time::Instant;
use std::{env, process};

use getopts::Options;

//-----------------------------------------------------------------------------

fn main() -> Result<(), String> {
    let start = Instant::now();
    let config = Config::new();

    eprintln!("Loading GBZ graph {}", config.gbz_filename);
    let graph: GBZ = serialize::load_from(&config.gbz_filename).map_err(|e| format!("Failed to load GBZ graph: {}", e))?;

    eprintln!("Building chains");
    let construction_start = Instant::now();
    let chains = algorithms::find_chains(&graph);
    let seconds = construction_start.elapsed().as_secs_f64();
    eprintln!("Built chains in {:.3} seconds", seconds);
    serialize::serialize_to(&chains, &config.output_filename).map_err(|e| format!("Failed to serialize chains: {}", e))?;
    eprintln!("There are {} chains with a total of {} links", chains.len(), chains.links());

    if let Some(truth_filename) = config.truth_filename.as_ref() {
        eprintln!("Loading truth chains from {}", truth_filename);
        let truth: Chains = serialize::load_from(truth_filename).map_err(|e| format!("Failed to load truth chains: {}", e))?;
        let mut identical = true;

        if chains.len() != truth.len() {
            eprintln!("Number of chains differ: {} vs {}", chains.len(), truth.len());
            identical = false;
        }
        if chains.links() != truth.links() {
            eprintln!("Number of links differ: {} vs {}", chains.links(), truth.links());
            identical = false;
        }

        let mut chains_iter = chains.iter().peekable();
        let mut truth_iter = truth.iter().peekable();
        while chains_iter.peek().is_some() || truth_iter.peek().is_some() {
            match (chains_iter.peek(), truth_iter.peek()) {
                (Some(chain), Some(truth)) => {
                    if chain == truth {
                        chains_iter.next();
                        truth_iter.next();
                    } else if chain < truth {
                        if support::encoded_edge_is_canonical(chain.0, chain.1) {
                            eprintln!(
                                "Link in found chains but not in truth: {:?}: ({} {}, {} {})",
                                chain, support::node_id(chain.0), support::node_orientation(chain.0),
                                support::node_id(chain.1), support::node_orientation(chain.1)
                            );
                            identical = false;
                        }
                        chains_iter.next();
                    } else {
                        if support::encoded_edge_is_canonical(truth.0, truth.1) {
                            eprintln!(
                                "Link in truth but not in found chains: {:?}: ({} {}, {} {})",
                                truth, support::node_id(truth.0), support::node_orientation(truth.0),
                                support::node_id(truth.1), support::node_orientation(truth.1)
                            );
                            identical = false;
                        }
                        truth_iter.next();
                    }
                }
                (Some(chain), None) => {
                    if support::encoded_edge_is_canonical(chain.0, chain.1) {
                        eprintln!(
                            "Link in found chains but not in truth: {:?}: ({} {}, {} {})",
                            chain, support::node_id(chain.0), support::node_orientation(chain.0),
                            support::node_id(chain.1), support::node_orientation(chain.1)
                        );
                        identical = false;
                    }
                    chains_iter.next();
                }
                (None, Some(truth)) => {
                    if support::encoded_edge_is_canonical(truth.0, truth.1) {
                        eprintln!(
                            "Link in truth but not in found chains: {:?}: ({} {}, {} {})",
                            truth, support::node_id(truth.0), support::node_orientation(truth.0),
                            support::node_id(truth.1), support::node_orientation(truth.1)
                        );
                        identical = false;
                    }
                    truth_iter.next();
                }
                (None, None) => unreachable!(),
            }
        }

        if identical {
            eprintln!("Built chains are identical to the truth");
        }
    }

    eprintln!("Used {:.3} seconds", start.elapsed().as_secs_f64());
    internal::report_memory_usage();
    eprintln!();
    Ok(())
}

//-----------------------------------------------------------------------------

pub struct Config {
    pub gbz_filename: String,
    pub output_filename: String,
    pub truth_filename: Option<String>,
}

impl Config {
    pub fn new() -> Config {
        let args: Vec<String> = env::args().collect();
        let program = args[0].clone();

        let mut opts = Options::new();
        opts.optflag("h", "help", "print this help");
        opts.optopt("t", "truth", "existing chains file that is assumed to be correct", "FILE");
        let matches = match opts.parse(&args[1..]) {
            Ok(m) => m,
            Err(f) => {
                eprintln!("{}", f);
                process::exit(1);
            }
        };

        let mut config = Config {
            gbz_filename: String::new(),
            output_filename: String::new(),
            truth_filename: None,
        };

        if matches.opt_present("h") {
            let header = format!("Usage: {} [options] graph.gbz graph.chains", program);
            eprint!("{}", opts.usage(&header));
            process::exit(0);
        }
        if let Some(arg) = matches.opt_str("t") {
            config.truth_filename = Some(arg);
        }

        if matches.free.len() != 2 {
            let header = format!("Usage: {} [options] graph.gbz graph.chains", program);
            eprint!("{}", opts.usage(&header));
            process::exit(1);
        }
        config.gbz_filename = matches.free[0].clone();
        config.output_filename = matches.free[1].clone();

        config
    }
}

//-----------------------------------------------------------------------------

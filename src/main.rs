mod ce;
mod structure;
mod visualization;
extern crate nalgebra as na;
use crate::structure::Geometry;
use clap::Parser;
use log::{debug, error, info, warn};
use std::path::Path;
use std::process::exit;
use std::time::Instant;
use structure::Validations;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short = 'm', help = "Path to the mobile structure")]
    mobile: String,

    #[arg(short = 't', help = "Path to the target structure")]
    target: String,

    #[arg(
        long = "output",
        help = "Save aligned structures as PDB files; they will be save in the same directory as `_aln.pdb`",
        default_value = "false"
    )]
    output: bool,

    #[arg(
        short = 'v',
        long = "verbose",
        help = "Increase output verbosity",
        default_value = "false"
    )]
    verbose: bool,

    #[arg(
        long = "randomize",
        help = "Randomly rotate both the mobile and the target for developlment/debug purposes",
        default_value = "false"
    )]
    randomize: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Initialize logger
    if args.verbose {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    }

    info!("Starting protein structure alignment");
    debug!("Verbose mode enabled");

    // Create a folder called `debug`, to save the intermediary files
    match std::fs::create_dir_all("debug") {
        Ok(_) => debug!("Created 'debug' directory"),
        Err(e) => warn!("Failed to create 'debug' directory: {}", e),
    }

    let mut pdb_i = match pdbtbx::open(&args.mobile, pdbtbx::StrictnessLevel::Medium) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            error!("Error opening PDB file {:?}: {:?}", args.mobile, e);
            exit(1)
        }
    };

    let mut pdb_j = match pdbtbx::open(&args.target, pdbtbx::StrictnessLevel::Medium) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            error!("Error opening PDB file {:?}: {:?}", args.target, e);
            exit(1)
        }
    };

    // Apply validations
    if pdb_i.is_multimodel() || pdb_j.is_multimodel() {
        error!("Multimodel PDB files are not supported");
        exit(1);
    }

    if args.randomize {
        warn!("Randomly rotating structures for development purposes");
        pdb_i.randomly_rotate();
        pdb_j.randomly_rotate();
    }

    // Take note of how long it took
    let start = Instant::now();
    let (pdb_i_aln, pdb_j_aln) = ce::align(pdb_i, pdb_j);
    let duration = start.elapsed();
    info!("Alignment complete, took: {:?}", duration);

    // If output option is set, save the aligned structures
    if args.output {
        let mobile_path = Path::new(&args.mobile);
        let target_path = Path::new(&args.target);

        let mobile_stem = mobile_path.file_stem().unwrap().to_str().unwrap();
        let target_stem = target_path.file_stem().unwrap().to_str().unwrap();

        let mobile_output = format!("{}_aln.pdb", mobile_stem);
        let target_output = format!("{}_aln.pdb", target_stem);

        let _ = pdbtbx::save_pdb(&pdb_i_aln, &mobile_output, pdbtbx::StrictnessLevel::Strict);
        let _ = pdbtbx::save_pdb(&pdb_j_aln, &target_output, pdbtbx::StrictnessLevel::Strict);

        info!(
            "Saved aligned structures as {} and {}",
            mobile_output, target_output
        );
    }

    info!("Alignment process completed successfully");

    Ok(())
}

mod ce;
mod structure;
mod utils;
mod visualization;
use clap::Parser;
use log::{debug, error, info};
use std::process::exit;
use std::time::Instant;
use structure::{Geometry, Validations};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short = 'm', help = "Path to the mobile structure")]
    mobile: String,

    #[arg(short = 't', help = "Path to the target structure")]
    target: String,

    #[arg(
        short = 'o',
        long = "output",
        help = "Save aligned structures as PDB files; they will be saved in the current directory as `<name>_aln.pdb`",
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
        short = 'r',
        long = "randomize",
        help = "Randomly rotate both the mobile and the target for development/debug purposes",
        default_value = "false"
    )]
    randomize: bool,

    #[arg(
        short = 'p',
        long = "plot",
        help = "Save the alignment path plot as plot.png",
        default_value = "false"
    )]
    plot: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.verbose {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
        debug!("Verbose mode enabled");
    } else {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    }

    let mut pdb_i = match pdbtbx::open(&args.mobile) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            error!("Error opening PDB file {:?}: {:?}", args.mobile, e);
            exit(1)
        }
    };

    let mut pdb_j = match pdbtbx::open(&args.target) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            error!("Error opening PDB file {:?}: {:?}", args.target, e);
            exit(1)
        }
    };

    if pdb_i.is_multimodel() || pdb_j.is_multimodel() {
        error!("Multimodel PDB files are not supported");
        exit(1);
    }

    if args.randomize {
        info!("Randomly rotating structures for development purposes");
        pdb_i.randomly_rotate();
        pdb_j.randomly_rotate();
    }

    let start = Instant::now();
    let (pdb_i_aln, pdb_j_aln, aligned_rmsd, _) = ce::align(pdb_i, pdb_j, args.plot);
    let _ = start.elapsed();

    println!("{:.3}", aligned_rmsd);

    if args.plot {
        println!("Saved: plot.png");
    }

    if args.output {
        let mobile_output = utils::aln_filename(&args.mobile);
        let target_output = utils::aln_filename(&args.target);

        let _ = pdbtbx::save_pdb(&pdb_i_aln, &mobile_output, pdbtbx::StrictnessLevel::Strict);
        let _ = pdbtbx::save_pdb(&pdb_j_aln, &target_output, pdbtbx::StrictnessLevel::Strict);

        println!("Saved: {} and {}", mobile_output, target_output);
    }

    Ok(())
}

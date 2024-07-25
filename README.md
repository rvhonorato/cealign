# cealign

A command-line application implementing the [Combinatorial Extension Algorithm](https://doi.org/10.1093/protein/11.9.739) (`cealign`) for protein alignment, written in Rust.

## Project Status

⚠️ **Under Active Development** - This project is currently in development and **not** yet ready for production use.

## Description

`cealign` is a command-line tool that uses the Combinatorial Extension (CE) algorithm to align two protein structures. The CE algorithm is known for its effectiveness in structural alignment of proteins, especially when dealing with proteins that have low sequence similarity but similar structural features.

> There is another implementation of this algorithm is available in the widely used molecular visualizer [PyMol](https://pymolwiki.org/index.php/Cealign)

## Installation

```bash
git clone https://github.com/rvhonorato/cealign && cd cealign
cargo install --path .
```

## Usage

```bash
$ cealign -m data/2oob_A.pdb -t data/2oob_B.pdb
[2024-07-25T14:07:24Z INFO  cealign] Starting protein structure alignment
[2024-07-25T14:07:24Z INFO  cealign::ce] Initial RMSD (full): 4.978
[2024-07-25T14:07:24Z INFO  cealign::ce] Initial RMSD (fragment): 4.931
[2024-07-25T14:07:24Z INFO  cealign::ce] Final RMSD (fragment): 5.794
[2024-07-25T14:07:24Z INFO  cealign::ce] Final RMSD (full): 3.967
[2024-07-25T14:07:24Z INFO  cealign] Alignment complete, took: 22.594718ms
[2024-07-25T14:07:24Z INFO  cealign] Alignment process completed successfully
```

```text
$ cealign -h
Usage: cealign [OPTIONS] -m <MOBILE> -t <TARGET>

Options:
  -m <MOBILE>      Path to the mobile structure
  -t <TARGET>      Path to the target structure
      --output     Save aligned structures as PDB files; they will be save in the same directory as `_aln.pdb`
  -v, --verbose    Increase output verbosity
      --randomize  Randomly rotate both the mobile and the target for developlment/debug purposes
  -h, --help       Print help
  -V, --version    Print version
```

## Features

_Coming Soon_

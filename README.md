# cealign

[![ci](https://github.com/rvhonorato/cealign/actions/workflows/ci.yml/badge.svg)](https://github.com/rvhonorato/cealign/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/7a1c7929f01e4b5aa4d9c22bf7c704ee)](https://app.codacy.com/gh/rvhonorato/cealign/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![crates.io](https://img.shields.io/crates/v/cealign.svg)](https://crates.io/crates/cealign)
[![docs.rs](https://docs.rs/cealign/badge.svg)](https://docs.rs/cealign)
[![DOI](https://zenodo.org/badge/816861096.svg)](https://doi.org/10.5281/zenodo.18892597)

A Rust implementation of the [Combinatorial Extension (CE) algorithm](https://doi.org/10.1093/protein/11.9.739) for protein structure alignment, available as both a command-line tool and a library.

## Why cealign?

CE performs **local structural alignment**: it builds a similarity matrix from all possible aligned fragment pairs (AFPs) between the two structures and finds the optimal path through it — without requiring a pre-defined residue correspondence. This is different from global superposition methods (e.g., Kabsch/RMSD fitting), which force every residue to be paired and therefore require structures of identical length and topology.

Because the alignment is local, CE handles three situations where global methods fail or produce misleading results:

1. **Different lengths** — you can align a 46-residue protein against a 153-residue one; only the structurally conserved core is superimposed.
2. **Low or no sequence similarity** — CE is particularly effective in the "twilight zone" (≲20% sequence identity) where sequence-based alignment methods are unreliable, because it operates purely on Cα geometry.
3. **Partial structural matches** — if two proteins share only a domain or a motif, CE will identify and align just that region.

The reported RMSD is computed exclusively over the `n_aligned` Cα atoms in the final superposition. Insertions, deletions, and divergent loops are excluded. A global RMSD over all residues would be undefined — or misleading — in the same scenarios.

## Usage

### CLI

```
Usage: cealign [OPTIONS] -m <MOBILE> -t <TARGET>

Options:
  -m <MOBILE>   Path to the mobile structure
  -t <TARGET>   Path to the target structure
  -o, --output  Save aligned structures as PDB files (`<name>_aln.pdb`)
  -v, --verbose Increase output verbosity
  -p, --plot    Save the alignment path plot as plot.png [requires --features plot]
  -h, --help    Print help
  -V, --version Print version
```

Align two structures and print the RMSD:

```bash
cealign -m mobile.pdb -t target.pdb
```

Save the aligned structures to disk:

```bash
cealign -m mobile.pdb -t target.pdb -o
```

Save the alignment path plot (requires `--features plot`):

```bash
cealign -m mobile.pdb -t target.pdb -p
```

### Library

```rust
use cealign::ce;

let (mobile, _) = pdbtbx::open("mobile.pdb").unwrap();
let (target, _) = pdbtbx::open("target.pdb").unwrap();

let (aligned_mobile, reference, rmsd, n_aligned) = ce::align(mobile, target, false);
println!("{:.3}", rmsd);
```

## Installation

### CLI

```bash
cargo install cealign
```

To include the optional alignment path plot feature:

```bash
cargo install cealign --features plot
```

> The `plot` feature requires `libfontconfig` on Linux (`sudo apt-get install libfontconfig-dev`).

### Library

Add to your `Cargo.toml`:

```toml
[dependencies]
cealign = "0.1"
```

## Examples

### Crambin vs thaumatin-like protein (X-ray vs NMR)

Align crambin X-ray structure ([1CRN](https://www.rcsb.org/structure/1CRN)) against the first NMR model of ([1CCM](https://www.rcsb.org/structure/1CCM)).

1CCM is a multi-model NMR ensemble, so we extract model 1 first:

```bash
# Download structures
curl -s https://files.rcsb.org/download/1CRN.pdb -o 1crn.pdb
curl -s https://files.rcsb.org/download/1CCM.pdb -o 1ccm_full.pdb

# Extract model 1 from the NMR ensemble
awk '/^MODEL/ && $2==1{p=1} p; /^ENDMDL/ && p{p=0}' 1ccm_full.pdb > 1ccm_1.pdb

# Align
cealign -m 1crn.pdb -t 1ccm_1.pdb
```

Expected output (RMSD in Å over the aligned Cα atoms):

```
1.331
```

Full run with verbose logging, alignment output, and plot (requires `--features plot`):

```bash
cealign -m 1crn.pdb -t 1ccm_1.pdb -v -o -p
```

```
[DEBUG cealign] Verbose mode enabled
[INFO  cealign::ce] Initial RMSD (full): 22.172
[DEBUG cealign::ce] find_path: 20 candidate paths in buffer
[DEBUG cealign::ce] Candidate path: 5 AFPs (40 residues), aligned RMSD = 1.331
...
[INFO  cealign::ce] Final RMSD (aligned CA, 40 residues): 1.331
[INFO  cealign::ce] Final RMSD (full): 1.921
1.331
Saved: plot.png
Saved: 1crn_aln.pdb and 1ccm_1_aln.pdb
```

The alignment path plot (`plot.png`) shows each AFP as a diagonal segment on a residue-index grid, mirroring Figure 2 from [Shindyalov & Bourne (1998)](https://doi.org/10.1093/protein/11.9.739):

![Alignment path plot](https://raw.githubusercontent.com/rvhonorato/cealign/refs/heads/main/assets/plot.png)

### Twilight-zone alignment: RSV integrase vs Mu transposase

This is the canonical example from the [PyMOL cealign wiki](https://pymolwiki.org/index.php/Cealign): aligning the catalytic core of Rous sarcoma virus integrase ([1C0M](https://www.rcsb.org/structure/1C0M), chain B) against bacteriophage Mu transposase ([1BCO](https://www.rcsb.org/structure/1BCO)). The two proteins share a conserved integrase-like fold despite very low sequence identity — a typical twilight-zone case where sequence-based alignment fails but structural alignment succeeds.

```bash
# Download structures
curl -s https://files.rcsb.org/download/1C0M.pdb -o 1c0m.pdb
curl -s https://files.rcsb.org/download/1BCO.pdb -o 1bco.pdb

# Extract chain B from 1C0M; keep only ATOM records (strip waters/ligands)
awk '/^ATOM/ && substr($0,22,1)=="B"' 1c0m.pdb > 1c0mB.pdb
grep "^ATOM" 1bco.pdb > 1bco_clean.pdb

# Align
cealign -m 1bco_clean.pdb -t 1c0mB.pdb -v
```

Expected output:

```
[INFO  cealign::ce] Initial RMSD (full): 96.092
[INFO  cealign::ce] Final RMSD (aligned CA, 152 residues): 4.958
[INFO  cealign::ce] Final RMSD (full): 25.534
4.958
```

152 residues are aligned between the two proteins despite their very low sequence identity — CE identified the shared integrase-like core and reports RMSD only over those structurally conserved residues.

| Superimposed | 1BCO (Mu transposase) | 1C0M chain B (RSV integrase) |
|:---:|:---:|:---:|
| ![Superimposed alignment](https://raw.githubusercontent.com/rvhonorato/cealign/refs/heads/main/assets/aln_1.png) | ![1BCO](https://raw.githubusercontent.com/rvhonorato/cealign/refs/heads/main/assets/aln_2.png) | ![1C0M chain B](https://raw.githubusercontent.com/rvhonorato/cealign/refs/heads/main/assets/aln_3.png) |

## Limitations

- Single-model PDB files only — multi-model NMR ensembles must be split before use
- Alignment is performed on Cα atoms; all heavy atoms are transformed using the resulting rotation and translation

## See Also

- [PyMOL cealign](https://pymolwiki.org/index.php/Cealign)
- [Biopython Bio.PDB.cealign](https://biopython.org/docs/dev/api/Bio.PDB.cealign.html)
- [Shindyalov & Bourne (1998)](https://doi.org/10.1093/protein/11.9.739) — original CE algorithm paper

## License

[0BSD](LICENSE)

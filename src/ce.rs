use crate::structure;
use crate::structure::Geometry;
use crate::visualization;
use log::{debug, info};
use std::collections::{HashMap, HashSet};

// Compute the inter-AFP distance D_ij (Equation 6, independent set) between
// an AFP already in the path at index `path_i`/`path_j` and a candidate AFP
// at index `j_a`/`j_b`. Uses the anti-diagonal pattern so each residue
// appears exactly once. Matches the inner loop in ccealignmodule.c findPath().
#[allow(clippy::too_many_arguments)]
fn inter_afp_distance(
    da: &HashMap<usize, HashMap<usize, f64>>,
    db: &HashMap<usize, HashMap<usize, f64>>,
    serials_a: &[usize],
    serials_b: &[usize],
    path_i: usize,
    path_j: usize,
    j_a: usize,
    j_b: usize,
    win_size: usize,
) -> f64 {
    let mut score = 0.0f64;
    // start-to-start
    score +=
        (da[&serials_a[path_i]][&serials_a[j_a]] - db[&serials_b[path_j]][&serials_b[j_b]]).abs();
    // end-to-end
    score += (da[&serials_a[path_i + win_size - 1]][&serials_a[j_a + win_size - 1]]
        - db[&serials_b[path_j + win_size - 1]][&serials_b[j_b + win_size - 1]])
        .abs();
    // anti-diagonal interior terms (k = 1..win_size-2)
    for k in 1..win_size - 1 {
        score += (da[&serials_a[path_i + k]][&serials_a[j_a + win_size - 1 - k]]
            - db[&serials_b[path_j + k]][&serials_b[j_b + win_size - 1 - k]])
            .abs();
    }
    score
}

// Returns up to MAX_KEPT candidate alignment paths as AFP index-pair lists.
// Mirrors the full findPath() logic in ccealignmodule.c including:
//   - winCache / allScoreBuffer / tIndex for the weighted total-score formula
//   - ring buffer of up to 20 best paths
//   - starting AFP requires S < D0; candidate AFPs require S <= D0
#[allow(unused_assignments)]
fn find_path(
    da: &HashMap<usize, HashMap<usize, f64>>,
    db: &HashMap<usize, HashMap<usize, f64>>,
    s: &HashMap<(usize, usize), f64>,
    win_size: usize,
) -> Vec<Vec<(usize, usize)>> {
    const D0: f64 = 3.0;
    const D1: f64 = 4.0;
    const MAX_KEPT: usize = 20;
    const GAP_MAX: usize = 30;

    let mut serials_a: Vec<usize> = da.keys().cloned().collect();
    serials_a.sort_unstable();
    let mut serials_b: Vec<usize> = db.keys().cloned().collect();
    serials_b.sort_unstable();

    let len_a = serials_a.len();
    let len_b = serials_b.len();
    let smaller = len_a.min(len_b);

    // Number of independent intra-AFP distance pairs per window
    let win_sum = (win_size - 1) * (win_size - 2) / 2;

    // winCache[i] = (i+1)*i*win_size/2 + (i+1)*win_sum
    // Cumulative weight for a path of length i+1 (matches C reference exactly).
    let win_cache: Vec<f64> = (0..=smaller)
        .map(|i| ((i + 1) * i * win_size / 2 + (i + 1) * win_sum) as f64)
        .collect();

    // Global best path (updated across all (iA, iB) starting points)
    let mut best_path_length = 0usize;
    let mut best_path_score = 1e6f64;
    let mut best_path: Vec<(usize, usize)> = vec![(usize::MAX, usize::MAX); smaller];

    // Ring buffer of up to MAX_KEPT paths (mirrors pathBuffer in C reference)
    let mut buffer_index = 0usize;
    let mut buffer_size = 0usize;
    let mut len_buffer = [0usize; MAX_KEPT];
    let mut score_buffer = [1e6f64; MAX_KEPT];
    let mut path_buffer: Vec<Vec<(usize, usize)>> = vec![Vec::new(); MAX_KEPT];

    // Shared state across all (iA, iB) starting points — mirrors allScoreBuffer
    // and tIndex in C. Each slot is always overwritten before it is read within
    // the same path traversal, so stale values from prior starting points are
    // harmless (same semantics as the C reference).
    let mut all_score_buffer = vec![vec![1e6f64; GAP_MAX * 2 + 1]; smaller + 1];
    let mut t_index = vec![0usize; smaller + 1];

    for i_a in 0..len_a {
        if best_path_length > 1 && i_a > len_a.saturating_sub(win_size * (best_path_length - 1)) {
            break;
        }

        for i_b in 0..len_b {
            let sa0 = serials_a[i_a];
            let sb0 = serials_b[i_b];

            // Starting AFP: strictly < D0  (C: `if (S[iA][iB] >= D0) continue`)
            let s0 = match s.get(&(sa0, sb0)) {
                Some(&v) if (0.0..D0).contains(&v) => v,
                _ => continue,
            };

            if best_path_length > 1 && i_b > len_b.saturating_sub(win_size * (best_path_length - 1))
            {
                break;
            }

            let mut cur_path = vec![(usize::MAX, usize::MAX); smaller];
            cur_path[0] = (i_a, i_b);
            let mut cur_path_length = 1usize;
            t_index[0] = 0;
            let mut cur_total_score = 0.0f64;

            loop {
                let (prev_i, prev_j) = cur_path[cur_path_length - 1];
                let mut gap_best_score = 1e6f64;
                let mut gap_best_index: Option<usize> = None;

                #[allow(clippy::needless_range_loop)]
                for g in 0..(GAP_MAX * 2 + 1) {
                    let j_a = prev_i + win_size + if (g + 1) % 2 == 0 { g.div_ceil(2) } else { 0 };
                    let j_b = prev_j + win_size + if (g + 1) % 2 != 0 { g.div_ceil(2) } else { 0 };

                    // Bounds: AFP at j_a needs atoms j_a..j_a+win_size-1
                    // (mirrors C: `if (jA > lenA-winSize-1)`)
                    if j_a > len_a.saturating_sub(win_size + 1)
                        || j_b > len_b.saturating_sub(win_size + 1)
                    {
                        continue;
                    }

                    let sa = serials_a[j_a];
                    let sb = serials_b[j_b];

                    // Candidate AFP: <= D0  (C: `if (S[jA][jB] > D0) continue`)
                    match s.get(&(sa, sb)) {
                        Some(&v) if (0.0..=D0).contains(&v) => {}
                        _ => continue,
                    }

                    // Inter-AFP score against every AFP already in the path
                    let mut cur_score = 0.0f64;
                    for &(pi, pj) in cur_path[..cur_path_length].iter() {
                        cur_score += inter_afp_distance(
                            da, db, &serials_a, &serials_b, pi, pj, j_a, j_b, win_size,
                        );
                    }
                    cur_score /= win_size as f64 * cur_path_length as f64;

                    if cur_score >= D1 {
                        continue;
                    }

                    if cur_score < gap_best_score {
                        cur_path[cur_path_length] = (j_a, j_b);
                        gap_best_score = cur_score;
                        gap_best_index = Some(g);
                        all_score_buffer[cur_path_length - 1][g] = cur_score;
                    }
                }

                match gap_best_index {
                    Some(g_best) => {
                        // Reconstruct (gA, gB) indices for the best candidate AFP.
                        let j_gap = g_best.div_ceil(2);
                        let (g_a, g_b) = if (g_best + 1) % 2 == 0 {
                            (prev_i + win_size + j_gap, prev_j + win_size)
                        } else {
                            (prev_i + win_size, prev_j + win_size + j_gap)
                        };

                        let s_ga_gb = s
                            .get(&(serials_a[g_a], serials_b[g_b]))
                            .copied()
                            .unwrap_or(0.0)
                            .max(0.0);

                        // score1: blend inter-AFP score with intra-AFP self-similarity
                        // of the new AFP (mirrors C: score1 formula).
                        let score1 = (all_score_buffer[cur_path_length - 1][g_best]
                            * win_size as f64
                            * cur_path_length as f64
                            + s_ga_gb * win_sum as f64)
                            / (win_size as f64 * cur_path_length as f64 + win_sum as f64);

                        // score2: combine score1 with running path total using
                        // winCache weights (mirrors C: score2 formula).
                        let prev_score = if cur_path_length > 1 {
                            all_score_buffer[cur_path_length - 2][t_index[cur_path_length - 1]]
                        } else {
                            s0
                        };

                        let score2 = (prev_score * win_cache[cur_path_length - 1]
                            + score1
                                * (win_cache[cur_path_length] - win_cache[cur_path_length - 1]))
                            / win_cache[cur_path_length];

                        cur_total_score = score2;

                        if cur_total_score > D1 {
                            break;
                        }

                        all_score_buffer[cur_path_length - 1][g_best] = cur_total_score;
                        t_index[cur_path_length] = g_best;
                        cur_path_length += 1;

                        // Update global best
                        if cur_path_length > best_path_length
                            || (cur_path_length == best_path_length
                                && cur_total_score < best_path_score)
                        {
                            best_path_length = cur_path_length;
                            best_path_score = cur_total_score;
                            best_path = cur_path.clone();
                        }
                    }
                    None => {
                        cur_path_length = cur_path_length.saturating_sub(1);
                        break;
                    }
                }
            }

            // Add global best to ring buffer if it improved vs. the current slot.
            // Mirrors the pathBuffer update logic in the C reference.
            if best_path_length > len_buffer[buffer_index]
                || (best_path_length == len_buffer[buffer_index]
                    && best_path_score < score_buffer[buffer_index])
            {
                buffer_index = if buffer_index == MAX_KEPT - 1 {
                    0
                } else {
                    buffer_index + 1
                };
                buffer_size = (buffer_size + 1).min(MAX_KEPT);

                let path_copy: Vec<(usize, usize)> = best_path[..best_path_length].to_vec();

                let store_idx = if buffer_index == 0 && buffer_size == MAX_KEPT {
                    MAX_KEPT - 1
                } else {
                    buffer_index - 1
                };

                path_buffer[store_idx] = path_copy;
                score_buffer[store_idx] = best_path_score;
                len_buffer[store_idx] = best_path_length;
            }
        }
    }

    if buffer_size == 0 {
        panic!("No alignment path found!");
    }

    debug!("find_path: {} candidate paths in buffer", buffer_size);
    path_buffer[..buffer_size].to_vec()
}

fn calculate_alignment_path(
    pdb_i: &pdbtbx::PDB,
    pdb_j: &pdbtbx::PDB,
    plot: bool,
) -> Vec<(Vec<usize>, Vec<usize>)> {
    let dm_i = structure::calc_distance_matrix(pdb_i);
    let dm_j = structure::calc_distance_matrix(pdb_j);

    let window_size = 8;
    let s = structure::calc_s(&dm_i, &dm_j, window_size);

    let mut serials_a: Vec<usize> = dm_i.keys().cloned().collect();
    serials_a.sort_unstable();
    let mut serials_b: Vec<usize> = dm_j.keys().cloned().collect();
    serials_b.sort_unstable();

    let candidate_paths = find_path(&dm_i, &dm_j, &s, window_size);

    // Try every candidate path; keep the expanded path with the best aligned
    // RMSD after Kabsch superposition — mirrors the PyMOL Python wrapper that
    // iterates pathBuffer and selects the lowest-RMSD solution.
    let mut best_expanded: Option<Vec<(Vec<usize>, Vec<usize>)>> = None;
    let mut best_index_path: Option<Vec<(usize, usize)>> = None;
    let mut best_rmsd = f64::MAX;

    for path in &candidate_paths {
        let serial_path: Vec<(Vec<usize>, Vec<usize>)> = path
            .iter()
            .map(|&(i, j)| {
                (
                    vec![serials_a[i], serials_a[i + window_size - 1]],
                    vec![serials_b[j], serials_b[j + window_size - 1]],
                )
            })
            .collect();

        let expanded = expand_path(serial_path, &dm_i, &dm_j);
        let (p, q) = structure::create_coordinate_vectors(&expanded, pdb_i, pdb_j);

        if p.len() != q.len() || p.is_empty() {
            continue;
        }

        let (rot, trans) = structure::calculate_transformation(&p, &q);
        let rmsd: f64 = {
            let sum_sq: f64 = p
                .iter()
                .zip(q.iter())
                .map(|(pi, qi)| (rot * pi + trans - qi).norm_squared())
                .sum();
            (sum_sq / p.len() as f64).sqrt()
        };

        debug!(
            "Candidate path: {} AFPs ({} residues), aligned RMSD = {:.3}",
            path.len(),
            p.len(),
            rmsd
        );

        if rmsd < best_rmsd {
            best_rmsd = rmsd;
            best_expanded = Some(expanded);
            best_index_path = Some(path.clone());
        }
    }

    let result = best_expanded.expect("No valid alignment path found!");
    result
        .iter()
        .for_each(|(a, b)| debug!("Best AFP: {:?} -> {:?}", a, b));

    // Generate alignment path plot (Figure 2 style from Shindyalov & Bourne 1998)
    if plot {
        if let Some(ref idx_path) = best_index_path {
            if let Err(e) =
                visualization::plot(idx_path, window_size, serials_a.len(), serials_b.len())
            {
                debug!("Could not write alignment plot: {}", e);
            }
        }
    }

    result
}

fn expand_path(
    path: Vec<(Vec<usize>, Vec<usize>)>,
    da: &HashMap<usize, HashMap<usize, f64>>,
    db: &HashMap<usize, HashMap<usize, f64>>,
) -> Vec<(Vec<usize>, Vec<usize>)> {
    let serials_a: HashSet<usize> = da.keys().cloned().collect();
    let serials_b: HashSet<usize> = db.keys().cloned().collect();

    path.into_iter()
        .map(|(range_a, range_b)| {
            let start_a = range_a[0].min(range_a[1]);
            let end_a = range_a[0].max(range_a[1]);
            let start_b = range_b[0].min(range_b[1]);
            let end_b = range_b[0].max(range_b[1]);

            let mut expanded_a: Vec<usize> = serials_a
                .iter()
                .filter(|&&x| x >= start_a && x <= end_a)
                .cloned()
                .collect();
            expanded_a.sort_unstable();

            let mut expanded_b: Vec<usize> = serials_b
                .iter()
                .filter(|&&x| x >= start_b && x <= end_b)
                .cloned()
                .collect();
            expanded_b.sort_unstable();

            (expanded_a, expanded_b)
        })
        .collect()
}

// Returns (aligned_mobile, reference, aligned_ca_rmsd, n_aligned_residues)
#[allow(non_snake_case)]
pub fn align(
    mut pdb_i: pdbtbx::PDB,
    mut pdb_j: pdbtbx::PDB,
    plot: bool,
) -> (pdbtbx::PDB, pdbtbx::PDB, f64, usize) {
    let rmsd = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Initial RMSD (full): {:.3}", rmsd);

    let pdb_i_all_atoms = pdb_i.clone();
    let pdb_j_all_atoms = pdb_j.clone();

    pdb_i.remove_atoms_by(|atom| atom.name() != "CA");
    pdb_j.remove_atoms_by(|atom| atom.name() != "CA");

    let path = calculate_alignment_path(&pdb_i, &pdb_j, plot);
    let n_aligned = path.iter().map(|(a, _)| a.len()).sum::<usize>();

    let (P, Q) = structure::create_coordinate_vectors(&path, &pdb_i, &pdb_j);

    let rmsd = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Initial RMSD (fragment): {:.3}", rmsd);

    let (rotation_matrix, translation_vector) = structure::calculate_transformation(&P, &Q);

    // RMSD over aligned residues only — matches PyMOL's reported metric.
    let aligned_ca_rmsd = {
        let sum_sq: f64 = P
            .iter()
            .zip(Q.iter())
            .map(|(p, q)| (rotation_matrix * p + translation_vector - q).norm_squared())
            .sum();
        (sum_sq / P.len() as f64).sqrt()
    };
    info!(
        "Final RMSD (aligned CA, {} residues): {:.3}",
        n_aligned, aligned_ca_rmsd
    );

    let mut pdb_i = pdb_i_all_atoms;
    let pdb_j = pdb_j_all_atoms;

    pdb_i.apply_rotation_and_translation(&rotation_matrix, &translation_vector);
    let final_rmsd_full = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Final RMSD (full): {:.3}", final_rmsd_full);

    (pdb_i, pdb_j, aligned_ca_rmsd, n_aligned)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Geometry;

    // Align 1CRN (crambin, X-ray) against each 1CCM NMR conformation.
    #[test]
    fn test_align_crambin_xray_vs_nmr() {
        // These base values were obtained by aligning `1crn` with `1ccm`
        //  using pymol's implementation of cealign
        let conformations = [
            ("data/1ccm_1.pdb", 1.331),
            ("data/1ccm_2.pdb", 1.156),
            ("data/1ccm_3.pdb", 0.921),
            ("data/1ccm_4.pdb", 1.193),
            ("data/1ccm_5.pdb", 1.067),
            ("data/1ccm_6.pdb", 1.027),
            ("data/1ccm_7.pdb", 1.066),
            ("data/1ccm_8.pdb", 1.034),
        ];

        for (path, pymol_rmsd) in conformations {
            let (reference, _) = pdbtbx::open("data/1crn.pdb", pdbtbx::StrictnessLevel::Medium)
                .expect("Failed to open data/1crn.pdb");

            let (mobile, _) = pdbtbx::open(path, pdbtbx::StrictnessLevel::Medium)
                .unwrap_or_else(|_| panic!("Failed to open {}", path));

            // Add a random rotation
            let mut mobile = mobile;
            mobile.randomly_rotate();

            let rmsd_before = structure::calc_rmsd(&reference, &mobile);
            let (_, _, rmsd_after, n_aligned) = align(mobile, reference, false);

            // println!(
            //     "{}: before: {:.3} Å  after (aligned CA, {} res): {:.3} Å  PyMOL: {:.3} Å",
            //     path, rmsd_before, n_aligned, rmsd_after, pymol_rmsd
            // );

            assert!(
                rmsd_after < rmsd_before,
                "{}: alignment should improve RMSD: before={:.3}, after={:.3}",
                path,
                rmsd_before,
                rmsd_after
            );

            assert_eq!(n_aligned, 40);
            assert!(
                rmsd_after < pymol_rmsd + 0.15,
                "{}: aligned CA RMSD {:.3} Å exceeds PyMOL {:.3} Å + 0.15 Å tolerance",
                path,
                rmsd_after,
                pymol_rmsd
            );
        }
    }
}

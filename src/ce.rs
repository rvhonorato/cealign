use crate::structure;
use crate::visualization;
use std::collections::HashMap;

use na::Vector3;

use crate::structure::Geometry;
use log::info;
use nalgebra::Matrix3;

fn find_path(
    da: &[Vec<f64>],
    db: &[Vec<f64>],
    s: &[Vec<f64>],
    win_size: usize,
) -> Vec<(Vec<usize>, Vec<usize>)> {
    // Note: This is the most important function! It's the core of the CE algorithm and can definitely be optimized and improved
    let mut path: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();

    let mut start_i = 0;
    let mut start_j = 0;
    let mut end_j = 0;

    for (i, _) in da.iter().enumerate() {
        if i < start_i {
            // This segment is already aligned, continue
            continue;
        }

        if i > da.len() - win_size {
            // we reached the end of the protein
            break;
        }

        let afp_a = vec![i, i + win_size];

        let mut possible_next_afp: HashMap<Vec<usize>, f64> = HashMap::new();

        for j in 0..db.len() {
            if j < start_j || j < end_j {
                // We have not walked far enough, continue
                continue;
            }

            if j > i + 30 {
                // We have walked too far, stop
                break;
            }

            let afp_b = vec![j, j + win_size];
            // println!("afp_a: {:?}, afp_b: {:?}", afp_a, afp_b);

            // How good is this alignment?
            let score = s[i][j];

            if score == -1.0 {
                // This alignment is not possible
                // println!(
                //     "No alignment possible for {:?} <-> {:?} - score = {}",
                //     afp_a, afp_b, score
                // );
                continue;
            }

            if score >= 3.0 {
                // This alignment is not good enough
                // println!(
                //     "- Discarded {:?} <-> {:?} - score = {:.2}",
                //     afp_a, afp_b, score
                // );
                continue;
            }

            // println!(
            //     "+ Accepted {:?} <-> {:?} - score = {:.2}",
            //     afp_a, afp_b, score
            // );

            // println!("afp_a: {:?}, afp_b: {:?}, score: {}", afp_a, afp_b, score);

            possible_next_afp.insert(afp_b, score);
        }

        // There might not be any possible next AFPs
        if possible_next_afp.is_empty() {
            // println!("## No possible next AFPs for {:?}", afp_a);
            continue;
        }

        // Which of the possible next AFPs is the best?
        let (best_next_afp, _score) = possible_next_afp
            .iter()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap();

        // println!(
        //     ">> Best: {:?} <-> {:?}, score: {:.2?}",
        //     afp_a, best_next_afp, score
        // );

        // Add the best next AFP to the path
        path.push((afp_a.clone(), best_next_afp.clone()));

        start_i = afp_a[1];
        start_j = best_next_afp[0] + 2;
        end_j = best_next_afp[1];

        // println!(
        //     "start_i: {}, start_j: {}, end_j: {}",
        //     start_i, start_j, end_j
        // );
    }

    // println!("{:?}", path);

    path
}

fn calculate_alignment_path(
    pdb_i: &pdbtbx::PDB,
    pdb_j: &pdbtbx::PDB,
) -> Vec<(Vec<usize>, Vec<usize>)> {
    // Calculate the distance matrixes, these are an all-vs-all distance matrix of the atoms of a given PDB
    let dm_i = structure::calc_distance_matrix(pdb_i);
    let dm_j = structure::calc_distance_matrix(pdb_j);

    // Calculate the similarity matrix
    let window_size = 8;
    let s = structure::calc_s(&dm_i, &dm_j, dm_i.len(), dm_j.len(), window_size);

    // Find the best path
    let path = find_path(&dm_i, &dm_j, &s, window_size);

    // Extra: visualize the path
    visualization::plot(path.clone(), dm_i.len(), dm_j.len()).unwrap();

    path
}

#[allow(non_snake_case)]
fn expand_path(
    path: &Vec<(Vec<usize>, Vec<usize>)>,
    pdb_i: &pdbtbx::PDB,
    pdb_j: &pdbtbx::PDB,
) -> (Vec<Vector3<f64>>, Vec<Vector3<f64>>) {
    // The `path` will tell us which pair of atoms need to be aligned with each other, use it to construct a vector of Points
    let mut P: Vec<Vector3<f64>> = Vec::new();
    let mut Q = Vec::new();

    for (p_a, p_b) in path {
        let expanded_i = (p_a[0]..p_a[1]).collect::<Vec<usize>>();
        let expanded_j = (p_b[0]..p_b[1]).collect::<Vec<usize>>();

        for (i, j) in expanded_i.iter().zip(expanded_j.iter()) {
            let atom_i = pdb_i.atoms().nth(*i).unwrap();
            let v = Vector3::new(atom_i.x(), atom_i.y(), atom_i.z());
            P.push(v);

            let atom_j = pdb_j.atoms().nth(*j).unwrap();
            let v = Vector3::new(atom_j.x(), atom_j.y(), atom_j.z());
            Q.push(v);
        }
    }

    (P, Q)
}

#[allow(non_snake_case)]
fn superimpose(
    mut pdb_i: pdbtbx::PDB,
    pdb_j: pdbtbx::PDB,
    reference_centroid: Vector3<f64>,
    U: Matrix3<f64>,
) -> (pdbtbx::PDB, pdbtbx::PDB) {
    for atom in pdb_i.atoms_mut() {
        let point = Vector3::new(atom.x(), atom.y(), atom.z());
        let centered_point = point - reference_centroid;
        let rotated_point = U * centered_point;
        let final_point = rotated_point + reference_centroid;

        let _ = atom.set_x(final_point.x);
        let _ = atom.set_y(final_point.y);
        let _ = atom.set_z(final_point.z);
    }

    // Apply final translation
    let centroid_i = pdb_i.calculate_centroid();
    let centroid_j = pdb_j.calculate_centroid();
    let translation = centroid_j - centroid_i;

    for atom in pdb_i.atoms_mut() {
        let _ = atom.set_x(atom.x() + translation.x);
        let _ = atom.set_y(atom.y() + translation.y);
        let _ = atom.set_z(atom.z() + translation.z);
    }

    // Calculate final RMSD
    let _ = pdbtbx::save_pdb(&pdb_i, "debug/mobile.pdb", pdbtbx::StrictnessLevel::Strict);
    let _ = pdbtbx::save_pdb(&pdb_j, "debug/target.pdb", pdbtbx::StrictnessLevel::Strict);
    let final_rmsd_full = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Final RMSD (full): {:.3}", final_rmsd_full);

    (pdb_i, pdb_j)
}

#[allow(non_snake_case)]
pub fn align(mut pdb_i: pdbtbx::PDB, mut pdb_j: pdbtbx::PDB) -> (pdbtbx::PDB, pdbtbx::PDB) {
    // Get the initial RMSD
    let rmsd = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Initial RMSD (full): {:.3}", rmsd);

    // Remove all atoms that are not `CA` from the PDB
    let pdb_i_all_atoms = pdb_i.clone();
    let pdb_j_all_atoms = pdb_j.clone();

    pdb_i.remove_atoms_by(|atom| atom.name() != "CA");
    pdb_j.remove_atoms_by(|atom| atom.name() != "CA");

    let _ = pdbtbx::save_pdb(&pdb_i, "debug/ca_i.pdb", pdbtbx::StrictnessLevel::Strict);
    let _ = pdbtbx::save_pdb(&pdb_j, "debug/ca_j.pdb", pdbtbx::StrictnessLevel::Strict);

    // Create the path using the CE algorithm
    let path = calculate_alignment_path(&pdb_i, &pdb_j);

    // Expand the path to get the points that need to be aligned
    let (P, Q) = expand_path(&path, &pdb_i, &pdb_j);

    let U = structure::calculate_rotation_matrix(&P, &Q);

    // Apply this rotation matrix to the CA only PDB
    let rmsd = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Initial RMSD (fragment): {:.3}", rmsd);
    pdb_i.apply_rotation(&U);

    let rmsd = structure::calc_rmsd(&pdb_i, &pdb_j);
    info!("Final RMSD (fragment): {:.3}", rmsd,);

    // Apply rotation to entire structure
    let pdb_i = pdb_i_all_atoms;
    let pdb_j = pdb_j_all_atoms;

    let centroid_P = structure::centroid(&P);
    superimpose(pdb_i, pdb_j, centroid_P, U)
}

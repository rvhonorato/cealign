use log::{debug, error};
use nalgebra::{Matrix3, Vector3};
use rand::Rng;
use std::collections::HashMap;
use std::process::exit;

pub trait Geometry {
    fn randomly_rotate(&mut self);
    fn apply_rotation_and_translation(
        &mut self,
        rotation_matrix: &Matrix3<f64>,
        translation_vector: &Vector3<f64>,
    );
}

impl Geometry for pdbtbx::PDB {
    fn apply_rotation_and_translation(
        &mut self,
        rotation_matrix: &Matrix3<f64>,
        translation_vector: &Vector3<f64>,
    ) {
        for atom in self.atoms_mut() {
            let point = Vector3::new(atom.x(), atom.y(), atom.z());

            // Apply rotation
            let rotated_point = rotation_matrix * point;

            // Apply translation
            let final_point = rotated_point + translation_vector;

            // Update atom coordinates
            let _ = atom.set_x(final_point.x);
            let _ = atom.set_y(final_point.y);
            let _ = atom.set_z(final_point.z);
        }
    }

    fn randomly_rotate(&mut self) {
        let mut rng = rand::thread_rng();

        // Generate a random U matrix
        let angle_x = rng.gen_range(0.0..=2.0 * std::f64::consts::PI);
        let angle_y = rng.gen_range(0.0..=2.0 * std::f64::consts::PI);
        let angle_z = rng.gen_range(0.0..=2.0 * std::f64::consts::PI);

        let rot_x = Matrix3::new(
            1.0,
            0.0,
            0.0, //
            0.0,
            angle_x.cos(),
            -angle_x.sin(), //
            0.0,
            angle_x.sin(),
            angle_x.cos(), //
        );

        let rot_y = Matrix3::new(
            angle_y.cos(),
            0.0,
            angle_y.sin(), //
            0.0,
            1.0,
            0.0, //
            -angle_y.sin(),
            0.0,
            angle_y.cos(), //
        );

        let rot_z = Matrix3::new(
            angle_z.cos(),
            -angle_z.sin(),
            0.0, //
            angle_z.sin(),
            angle_z.cos(),
            0.0, //
            0.0,
            0.0,
            1.0, //
        );

        let rotation_matrix = rot_x * rot_y * rot_z;

        // Apply rotation and translations
        for atom in self.atoms_mut() {
            let position = Vector3::new(atom.x(), atom.y(), atom.z());

            // Apply rotation
            let rotated_position = rotation_matrix * position;

            // Update the coordinates
            let _ = atom.set_x(rotated_position.x);
            let _ = atom.set_y(rotated_position.y);
            let _ = atom.set_z(rotated_position.z);
        }
    }
}

#[allow(non_snake_case)]
pub fn calculate_transformation(
    P: &[Vector3<f64>],
    Q: &[Vector3<f64>],
) -> (Matrix3<f64>, Vector3<f64>) {
    let rotation_matrix = calculate_rotation_matrix(P, Q);

    let centroid_P = centroid(P);
    let centroid_Q = centroid(Q);

    let translation_vector = centroid_Q - rotation_matrix * centroid_P;

    (rotation_matrix, translation_vector)
}

pub trait Validations {
    fn is_multimodel(&self) -> bool;
}

impl Validations for pdbtbx::PDB {
    fn is_multimodel(&self) -> bool {
        self.models().count() > 1
    }
}

pub fn calc_distance_matrix(pdb: &pdbtbx::PDB) -> HashMap<usize, HashMap<usize, f64>> {
    pdb.atoms()
        .map(|atom1| {
            let serial1 = atom1.serial_number();
            (
                serial1,
                pdb.atoms()
                    .map(|atom2| {
                        let serial2 = atom2.serial_number();
                        (serial2, atom1.distance(atom2))
                    })
                    .collect::<HashMap<usize, f64>>(),
            )
        })
        .collect()
}

#[allow(non_snake_case)]
pub fn create_coordinate_vectors(
    expanded_path: &Vec<(Vec<usize>, Vec<usize>)>,
    pdb_a: &pdbtbx::PDB,
    pdb_b: &pdbtbx::PDB,
) -> (Vec<Vector3<f64>>, Vec<Vector3<f64>>) {
    let mut P: Vec<Vector3<f64>> = Vec::new();
    let mut Q: Vec<Vector3<f64>> = Vec::new();

    for (serials_a, serials_b) in expanded_path {
        for &serial_a in serials_a {
            if let Some(atom_a) = pdb_a.atoms().find(|atom| atom.serial_number() == serial_a) {
                let v = Vector3::new(atom_a.x(), atom_a.y(), atom_a.z());
                P.push(v);
            }
        }

        for &serial_b in serials_b {
            if let Some(atom_b) = pdb_b.atoms().find(|atom| atom.serial_number() == serial_b) {
                let v = Vector3::new(atom_b.x(), atom_b.y(), atom_b.z());
                Q.push(v);
            }
        }
    }

    (P, Q)
}

pub fn calc_s(
    d1: &HashMap<usize, HashMap<usize, f64>>,
    d2: &HashMap<usize, HashMap<usize, f64>>,
    win_size: usize,
) -> HashMap<(usize, usize), f64> {
    let mut serials_a: Vec<usize> = d1.keys().cloned().collect();
    serials_a.sort_unstable();
    let mut serials_b: Vec<usize> = d2.keys().cloned().collect();
    serials_b.sort_unstable();

    let len_a = serials_a.len();
    let len_b = serials_b.len();

    // Number of distance pairs per window: all (row, col) with col >= row + 2.
    // Matches the normalization from the original Vec-based implementation.
    let sum_size = ((win_size - 1) * (win_size - 2)) as f64 / 2.0;

    let mut scores = HashMap::new();

    for i_a in 0..len_a {
        if i_a + win_size > len_a {
            break;
        }
        for i_b in 0..len_b {
            if i_b + win_size > len_b {
                break;
            }

            let mut score = 0.0;

            for row in 0..win_size - 2 {
                for col in row + 2..win_size {
                    let sa_row = serials_a[i_a + row];
                    let sa_col = serials_a[i_a + col];
                    let sb_row = serials_b[i_b + row];
                    let sb_col = serials_b[i_b + col];

                    let dist_a = d1[&sa_row][&sa_col];
                    let dist_b = d2[&sb_row][&sb_col];

                    score += (dist_a - dist_b).abs();
                }
            }

            scores.insert((serials_a[i_a], serials_b[i_b]), score / sum_size);
        }
    }

    scores
}

pub fn calc_rmsd(receptor: &pdbtbx::PDB, ligand: &pdbtbx::PDB) -> f64 {
    let n = receptor.atoms().count();
    let sum_sq = ligand
        .atoms()
        .zip(receptor.atoms())
        .map(|(a, b)| a.distance(b).powi(2))
        .sum::<f64>();
    (sum_sq / n as f64).sqrt()
}

pub fn centroid(points: &[Vector3<f64>]) -> Vector3<f64> {
    points.iter().sum::<Vector3<f64>>() / points.len() as f64
}

#[allow(non_snake_case)]
pub fn kabsch(P: &[Vector3<f64>], Q: &[Vector3<f64>]) -> Matrix3<f64> {
    if P.len() != Q.len() || P.is_empty() {
        error!("The input sets are not of the same size or are empty.");
        debug!("P: {:?}", P);
        debug!("Q: {:?}", Q);
        exit(1)
    }

    let covariance_matrix = P
        .iter()
        .zip(Q.iter())
        .fold(Matrix3::zeros(), |acc, (p, q)| acc + p * q.transpose());

    let svd = covariance_matrix.svd(true, true);
    // TODO: Handle these unwraps
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    // Calculate the determinant of the product of U and V^T
    let det = (u * v_t).determinant();

    // Create the diagonal matrix
    let mut d = Matrix3::identity();
    if det < 0.0 {
        *d.index_mut((2, 2)) = -1.0;
    }

    // Calculate the rotation matrix
    let r = v_t.transpose() * d * u.transpose();

    // Ensure the determinant is positive (proper rotation)
    if r.determinant() < 0.0 {
        -r // Negate the entire matrix if it's a reflection
    } else {
        r
    }
}

#[allow(non_snake_case)]
pub fn calculate_rotation_matrix(P: &[Vector3<f64>], Q: &[Vector3<f64>]) -> Matrix3<f64> {
    let centroid_P = centroid(P);
    let centroid_Q = centroid(Q);

    // Move P and Q to their centers
    let P_centered: Vec<Vector3<f64>> = P.iter().map(|p| p - centroid_P).collect();
    let Q_centered: Vec<Vector3<f64>> = Q.iter().map(|q| q - centroid_Q).collect();

    kabsch(&P_centered, &Q_centered)
}

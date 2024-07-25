use rand::Rng;

extern crate nalgebra as na;
use log::error;
use na::{Matrix3, Vector3};
use std::process::exit;

pub trait Geometry {
    fn apply_rotation(&mut self, rotation_matrix: &Matrix3<f64>);
    fn randomly_rotate(&mut self);
    fn calculate_centroid(&self) -> Vector3<f64>;
}

impl Geometry for pdbtbx::PDB {
    fn calculate_centroid(&self) -> Vector3<f64> {
        let mut sum = Vector3::zeros();
        let mut count = 0;

        for atom in self.atoms() {
            sum += Vector3::new(atom.x(), atom.y(), atom.z());
            count += 1;
        }

        if count > 0 {
            sum / count as f64
        } else {
            Vector3::zeros()
        }
    }

    fn apply_rotation(&mut self, rotation_matrix: &Matrix3<f64>) {
        // Calculate the centroid
        // let centroid = self.calculate_centroid();

        // Apply rotation around the centroid
        for atom in self.atoms_mut() {
            let point = Vector3::new(atom.x(), atom.y(), atom.z());

            // Move to origin
            // let centered_point = point - centroid;

            // Apply rotation
            // let rotated_point = rotation_matrix * centered_point;
            let final_point = rotation_matrix * point;

            // Move back from origin
            // let final_point = rotated_point + centroid;

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

pub trait Validations {
    fn is_multimodel(&self) -> bool;
}

impl Validations for pdbtbx::PDB {
    fn is_multimodel(&self) -> bool {
        self.models().count() > 1
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

// Calculate an all-vs-all distance
pub fn calc_distance_matrix(pdb: &pdbtbx::PDB) -> Vec<Vec<f64>> {
    pdb.atoms()
        .map(|atom1| {
            pdb.atoms()
                .map(|atom2| atom1.distance(atom2))
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>()
}

pub fn calc_s(
    d1: &[Vec<f64>],
    d2: &[Vec<f64>],
    len_a: usize,
    len_b: usize,
    win_size: usize,
) -> Vec<Vec<f64>> {
    let sum_size = ((win_size - 1) * (win_size - 2)) as f64 / 2.0;

    let mut s = vec![vec![-1.0; len_b]; len_a];

    for i_a in 0..len_a {
        for i_b in 0..len_b {
            if i_a > len_a - win_size || i_b > len_b - win_size {
                continue;
            }

            let mut score = 0.0;

            for row in 0..win_size - 2 {
                for col in row + 2..win_size {
                    score += (d1[i_a + row][i_a + col] - d2[i_b + row][i_b + col]).abs();
                }
            }

            s[i_a][i_b] = score / sum_size;
        }
    }

    s
}

pub fn calc_rmsd(receptor: &pdbtbx::PDB, ligand: &pdbtbx::PDB) -> f64 {
    // The receptor does not move, so we don't need to calculate over it,
    // then the RMSD here is the L-RMSD

    // Calculate the sum of squared distances between corresponding atoms
    let n = receptor.atoms().count();
    let sum_squared_distances = ligand
        .atoms()
        .zip(receptor.atoms())
        .map(|(atom1, atom2)| atom1.distance(atom2))
        .sum::<f64>();
    (sum_squared_distances / n as f64).sqrt()
}

pub fn centroid(points: &[Vector3<f64>]) -> Vector3<f64> {
    points.iter().sum::<Vector3<f64>>() / points.len() as f64
}

#[allow(non_snake_case)]
pub fn kabsch(P: &[Vector3<f64>], Q: &[Vector3<f64>]) -> Matrix3<f64> {
    if P.len() != Q.len() || P.is_empty() {
        error!("The input sets are not of the same size or are empty.");
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

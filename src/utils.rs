use std::path::Path;

/// Derive the output filename for an aligned structure.
///
/// Strips the directory and extension from `path` and appends `_aln.pdb`.
///
/// # Examples
///
/// ```
/// assert_eq!(cealign::utils::aln_filename("structures/mobile.pdb"), "mobile_aln.pdb");
/// ```
pub fn aln_filename(path: &str) -> String {
    let stem = Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    format!("{}_aln.pdb", stem)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aln_filename_simple() {
        assert_eq!(aln_filename("structure.pdb"), "structure_aln.pdb");
    }

    #[test]
    fn test_aln_filename_with_path() {
        assert_eq!(aln_filename("/some/dir/mobile.pdb"), "mobile_aln.pdb");
    }

    #[test]
    fn test_aln_filename_no_extension() {
        assert_eq!(aln_filename("structure"), "structure_aln.pdb");
    }
}

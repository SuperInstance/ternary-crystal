//! # ternary-crystal
//! 
//! Crystallography and lattice symmetry in ternary space.
//! Point groups, space groups, Brillouin zones, and diffraction patterns.

#![forbid(unsafe_code)]
#![no_std]

extern crate alloc;
use alloc::vec::Vec;
use alloc::vec;

/// A 2D or 3D lattice point in ternary coordinates
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LatticePoint {
    pub coords: [i8; 3],
}

impl LatticePoint {
    pub fn new(x: i8, y: i8, z: i8) -> Self {
        Self { coords: [((x % 3) + 3) % 3, ((y % 3) + 3) % 3, ((z % 3) + 3) % 3] }
    }

    /// Dot product modulo 3
    pub fn dot(&self, other: &Self) -> i8 {
        let raw: i8 = self.coords.iter()
            .zip(other.coords.iter())
            .map(|(a, b)| a * b)
            .sum();
        ((raw % 3) + 3) % 3
    }

    /// Distance in ternary metric (sum of absolute differences mod 3)
    pub fn ternary_distance(&self, other: &Self) -> i8 {
        self.coords.iter()
            .zip(other.coords.iter())
            .map(|(a, b)| {
                let d = (a - b).abs();
                d.min(3 - d) // ternary distance wraps
            })
            .sum()
    }
}

/// Symmetry operation on ternary lattice
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetryOp {
    Identity,
    Inversion,
    /// Rotation by 120 degrees around axis (0=x, 1=y, 2=z)
    C3(usize),
    /// Mirror plane perpendicular to axis
    Mirror(usize),
}

impl SymmetryOp {
    /// Apply symmetry operation to a lattice point
    pub fn apply(&self, point: &LatticePoint) -> LatticePoint {
        let [x, y, z] = point.coords;
        match self {
            SymmetryOp::Identity => LatticePoint::new(x, y, z),
            SymmetryOp::Inversion => LatticePoint::new(-x, -y, -z),
            SymmetryOp::C3(0) => LatticePoint::new(x, z, -y), // rotate around x
            SymmetryOp::C3(1) => LatticePoint::new(-z, y, x), // rotate around y
            SymmetryOp::C3(2) => LatticePoint::new(y, -x, z),
            SymmetryOp::C3(_) => *point,
            SymmetryOp::Mirror(0) => LatticePoint::new(-x, y, z),
            SymmetryOp::Mirror(1) => LatticePoint::new(x, -y, z),
            SymmetryOp::Mirror(2) => LatticePoint::new(x, y, -z),
            SymmetryOp::Mirror(_) => *point,
        }
    }

    /// Check if a lattice point is invariant under this operation
    pub fn fixes(&self, point: &LatticePoint) -> bool {
        &self.apply(point) == point
    }
}

/// A point group — the set of symmetry operations that leave a point fixed
#[derive(Debug, Clone)]
pub struct PointGroup {
    pub operations: Vec<SymmetryOp>,
    pub name: &'static str,
}

impl PointGroup {
    /// C1: trivial group (identity only)
    pub fn c1() -> Self {
        Self { operations: vec![SymmetryOp::Identity], name: "C1" }
    }

    /// C3i: inversion + C3 rotations
    pub fn c3i() -> Self {
        Self {
            operations: vec![
                SymmetryOp::Identity,
                SymmetryOp::Inversion,
                SymmetryOp::C3(2),
            ],
            name: "C3i",
        }
    }

    /// Full cubic group for ternary lattice
    pub fn ternary_cubic() -> Self {
        Self {
            operations: vec![
                SymmetryOp::Identity,
                SymmetryOp::Inversion,
                SymmetryOp::C3(0),
                SymmetryOp::C3(1),
                SymmetryOp::C3(2),
                SymmetryOp::Mirror(0),
                SymmetryOp::Mirror(1),
                SymmetryOp::Mirror(2),
            ],
            name: "Oh",
        }
    }

    /// Count fixed points of a lattice under this group
    pub fn fixed_point_count(&self, lattice_size: usize) -> usize {
        let mut count = 0;
        for x in 0..lattice_size {
            for y in 0..lattice_size {
                for z in 0..lattice_size {
                    let p = LatticePoint::new(x as i8, y as i8, z as i8);
                    if self.operations.iter().all(|op| op.fixes(&p)) {
                        count += 1;
                    }
                }
            }
        }
        count
    }
}

/// A crystal unit cell
#[derive(Debug, Clone)]
pub struct UnitCell {
    pub basis: Vec<LatticePoint>,
    pub lattice_vectors: [LatticePoint; 3],
}

impl UnitCell {
    /// Simple cubic ternary unit cell
    pub fn simple_cubic() -> Self {
        Self {
            basis: vec![LatticePoint::new(0, 0, 0)],
            lattice_vectors: [
                LatticePoint::new(1, 0, 0),
                LatticePoint::new(0, 1, 0),
                LatticePoint::new(0, 0, 1),
            ],
        }
    }

    /// Body-centered ternary unit cell
    pub fn body_centered() -> Self {
        Self {
            basis: vec![
                LatticePoint::new(0, 0, 0),
                LatticePoint::new(1, 1, 1),
            ],
            lattice_vectors: [
                LatticePoint::new(1, 0, 0),
                LatticePoint::new(0, 1, 0),
                LatticePoint::new(0, 0, 1),
            ],
        }
    }

    /// Generate all lattice points within a radius
    pub fn generate(&self, radius: i8) -> Vec<LatticePoint> {
        let mut points = Vec::new();
        for i in -radius..=radius {
            for j in -radius..=radius {
                for k in -radius..=radius {
                    for b in &self.basis {
                        let px = (i * 1 + b.coords[0] + 3) % 3;
                        let py = (j * 1 + b.coords[1] + 3) % 3;
                        let pz = (k * 1 + b.coords[2] + 3) % 3;
                        points.push(LatticePoint::new(px, py, pz));
                    }
                }
            }
        }
        points
    }
}

/// Brillouin zone — the Wigner-Seitz cell in reciprocal space
pub fn brillouin_zone(lattice: &[LatticePoint; 3]) -> Vec<LatticePoint> {
    // For ternary, the Brillouin zone is the set of reciprocal lattice points
    // within the first Brillouin zone
    let mut zone = Vec::new();
    for x in 0i8..3 {
        for y in 0i8..3 {
            for z in 0i8..3 {
                let p = LatticePoint::new(x, y, z);
                // A point is in the BZ if it's closer to origin than any other lattice point
                let d_origin = p.ternary_distance(&LatticePoint::new(0, 0, 0));
                let mut in_zone = true;
                for lv in lattice {
                    let shifted = LatticePoint::new(
                        (p.coords[0] - lv.coords[0] + 3) % 3,
                        (p.coords[1] - lv.coords[1] + 3) % 3,
                        (p.coords[2] - lv.coords[2] + 3) % 3,
                    );
                    if shifted.ternary_distance(&LatticePoint::new(0, 0, 0)) < d_origin {
                        in_zone = false;
                        break;
                    }
                }
                if in_zone {
                    zone.push(p);
                }
            }
        }
    }
    zone
}

/// Diffraction pattern — structure factor for a ternary crystal
pub fn structure_factor(cell: &UnitCell, q: &LatticePoint) -> i8 {
    // F(q) = sum over basis of exp(2*pi*i*q·r) mod 3
    // In ternary: exp(2*pi*i * dot(q,r)/3) maps to powers of the 3rd root of unity
    let mut phase_sum: i8 = 0;
    for atom in &cell.basis {
        let dot = q.dot(atom) % 3;
        phase_sum = (phase_sum + dot) % 3;
    }
    phase_sum
}

/// Count distinct crystal structures modulo symmetry
pub fn count_orbits(group: &PointGroup, size: usize) -> usize {
    // Burnside's lemma: |orbits| = (1/|G|) * sum_{g in G} |Fix(g)|
    let mut total_fixed = 0usize;
    for op in &group.operations {
        let mut fixed = 0usize;
        for x in 0..size {
            for y in 0..size {
                for z in 0..size {
                    let p = LatticePoint::new(x as i8, y as i8, z as i8);
                    if op.fixes(&p) {
                        fixed += 1;
                    }
                }
            }
        }
        total_fixed += fixed;
    }
    total_fixed / group.operations.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lattice_point_creation() {
        let p = LatticePoint::new(1, -1, 2);
        assert_eq!(p.coords, [1, 2, 2]); // -1 mod 3 = 2
    }

    #[test]
    fn test_dot_product() {
        let a = LatticePoint::new(1, 1, 0);
        let b = LatticePoint::new(1, 0, 1);
        assert_eq!(a.dot(&b), 1);
    }

    #[test]
    fn test_ternary_distance() {
        let a = LatticePoint::new(0, 0, 0);
        let b = LatticePoint::new(1, 1, 1);
        assert_eq!(a.ternary_distance(&b), 3);
        // Distance wraps: -1 ≡ 2, so distance from 0 to 2 is min(2, 1) = 1
        let c = LatticePoint::new(2, 0, 0);
        assert_eq!(a.ternary_distance(&c), 1);
    }

    #[test]
    fn test_symmetry_identity() {
        let p = LatticePoint::new(1, -1, 0);
        let result = SymmetryOp::Identity.apply(&p);
        assert_eq!(result.coords, [1, 2, 0]);
    }

    #[test]
    fn test_symmetry_inversion() {
        let p = LatticePoint::new(1, 0, -1);
        let result = SymmetryOp::Inversion.apply(&p);
        assert_eq!(result.coords, [2, 0, 1]);
    }

    #[test]
    fn test_c3_rotation_z() {
        // C3 around z: (x,y,z) -> (y,-x,z)
        let p = LatticePoint::new(1, 0, 0);
        let result = SymmetryOp::C3(2).apply(&p);
        assert_eq!(result.coords[2], 0); // z unchanged
    }

    #[test]
    fn test_mirror_x() {
        let p = LatticePoint::new(1, 1, 1);
        let result = SymmetryOp::Mirror(0).apply(&p);
        assert_eq!(result.coords[0], 2); // x inverted
        assert_eq!(result.coords[1], 1); // y unchanged
    }

    #[test]
    fn test_point_group_c1_fixed_points() {
        let group = PointGroup::c1();
        // All points are fixed by identity
        assert_eq!(group.fixed_point_count(3), 27);
    }

    #[test]
    fn test_simple_cubic_unit_cell() {
        let cell = UnitCell::simple_cubic();
        assert_eq!(cell.basis.len(), 1);
    }

    #[test]
    fn test_body_centered_unit_cell() {
        let cell = UnitCell::body_centered();
        assert_eq!(cell.basis.len(), 2);
    }

    #[test]
    fn test_generate_lattice() {
        let cell = UnitCell::simple_cubic();
        let points = cell.generate(1);
        assert!(!points.is_empty());
    }

    #[test]
    fn test_structure_factor() {
        let cell = UnitCell::simple_cubic();
        let q = LatticePoint::new(0, 0, 0);
        // F(0) = sum of phases = 0 (single atom at origin, dot = 0)
        assert_eq!(structure_factor(&cell, &q), 0);
    }

    #[test]
    fn test_structure_factor_bcc() {
        let cell = UnitCell::body_centered();
        let q = LatticePoint::new(1, 1, 1);
        // For BCC: F(111) = 1 + exp(2πi*3/3) = 1 + 1 = 2 in ternary → 2
        let sf = structure_factor(&cell, &q);
        // q·r1 = 0, q·r2 = 1+1+1 = 3 ≡ 0 mod 3. Sum = 0
        assert_eq!(sf, 0);
    }

    #[test]
    fn test_count_orbits() {
        let group = PointGroup::c1();
        let orbits = count_orbits(&group, 2);
        // With only identity, every point is its own orbit
        assert_eq!(orbits, 8); // 2^3
    }

    #[test]
    fn test_brillouin_zone() {
        let lattice = [
            LatticePoint::new(1, 0, 0),
            LatticePoint::new(0, 1, 0),
            LatticePoint::new(0, 0, 1),
        ];
        let zone = brillouin_zone(&lattice);
        assert!(!zone.is_empty());
    }
}

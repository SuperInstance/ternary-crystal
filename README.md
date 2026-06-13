# Ternary Crystal

Crystallography and **lattice symmetry in ternary space** — point groups, space groups, Brillouin zones, and diffraction patterns for ternary {-1, 0, +1} lattices. All operations use modular arithmetic modulo 3 (GF(3)), making the ternary lattice a natural object for group-theoretic analysis.

## Why It Matters

Classical crystallography uses real-valued coordinates with 230 space groups and 32 point groups. A ternary lattice — where every coordinate is in {0, 1, 2} (mod 3) — is a finite geometry with exactly $3^3 = 27$ points per unit cell. This radically constrains the symmetry landscape:

- **Finite group orders**: Every symmetry group acts on at most 27 points
- **GF(3) arithmetic**: Dot products, distances, and structure factors use modular arithmetic
- **Clean enumeration**: Burnside's lemma gives exact orbit counts

The ternary lattice is not just a mathematical curiosity — it models discrete structures where ternary {-1, 0, +1} values represent material properties (e.g., spin states, charge configurations, compositional ordering in high-entropy alloys).

## How It Works

### Lattice Points

Every point lives in $\mathbb{Z}_3^3$ — coordinates modulo 3:

$$(x, y, z) \in (\mathbb{Z}/3\mathbb{Z})^3$$

Construction normalizes: `LatticePoint::new(x, y, z)` applies `((x % 3) + 3) % 3` to each coordinate.

### Distance Metric

The ternary distance between points is:

$$d(\mathbf{a}, \mathbf{b}) = \sum_{i=1}^{3} \min(|a_i - b_i|, \; 3 - |a_i - b_i|)$$

This wraps around: the distance from 0 to 2 is $\min(2, 1) = 1$, not 2 — because in a ternary ring, 2 ≡ -1.

### Dot Product (mod 3)

$$\mathbf{a} \cdot \mathbf{b} = \left(\sum_{i=1}^{3} a_i b_i\right) \bmod 3$$

### Symmetry Operations

| Operation | Action on (x, y, z) |
|-----------|---------------------|
| Identity | (x, y, z) |
| Inversion | (-x, -y, -z) |
| C₃ rotation (axis z) | (y, -x, z) |
| Mirror (⊥ to x) | (-x, y, z) |

A point is **fixed** by an operation if applying the operation yields the same point (mod 3).

### Burnside's Lemma

The number of distinct crystal structures (orbits) under a group $G$ acting on the lattice is:

$$|\text{Orbits}| = \frac{1}{|G|} \sum_{g \in G} |\text{Fix}(g)|$$

where $|\text{Fix}(g)|$ is the number of lattice points fixed by operation $g$. For the trivial group C₁ acting on a size-2 lattice: $|G| = 1$, $|\text{Fix}(e)| = 2^3 = 8$ orbits.

### Structure Factor

For diffraction from a unit cell with basis $\{r_j\}$ at scattering vector $\mathbf{q}$:

$$F(\mathbf{q}) = \sum_{j} e^{2\pi i \mathbf{q} \cdot \mathbf{r}_j / 3} \pmod{3}$$

In GF(3), the exponential maps to powers of the primitive 3rd root of unity $\omega = e^{2\pi i/3}$. The structure factor is the phase sum modulo 3.

### Complexity

| Operation | Time |
|-----------|------|
| `LatticePoint::dot(&other)` | O(1) |
| `SymmetryOp::apply(&point)` | O(1) |
| `PointGroup::fixed_point_count(n)` | O(G · n³) |
| `brillouin_zone(&lattice)` | O(27 · 3) = O(1) |
| `structure_factor(&cell, &q)` | O(B) — B = basis atoms |
| `count_orbits(&group, n)` | O(G · n³) |

Where G = group order, n = lattice linear size.

## Quick Start

```rust
use ternary_crystal::{LatticePoint, SymmetryOp, PointGroup, UnitCell, brillouin_zone, structure_factor, count_orbits};

// Create lattice points (coordinates are mod 3)
let p = LatticePoint::new(1, -1, 2);
assert_eq!(p.coords, [1, 2, 2]); // -1 mod 3 = 2

// Apply symmetry
let inverted = SymmetryOp::Inversion.apply(&p);
let rotated  = SymmetryOp::C3(2).apply(&p); // 120° around z

// Point groups
let group = PointGroup::ternary_cubic(); // Full cubic group
let fixed = group.fixed_point_count(3);  // Points fixed by ALL operations

// Unit cells
let sc = UnitCell::simple_cubic();   // 1 basis atom
let bcc = UnitCell::body_centered();  // 2 basis atoms

// Diffraction
let q = LatticePoint::new(1, 1, 1);
let sf = structure_factor(&bcc, &q); // Structure factor at q

// Brillouin zone
let lattice = [
    LatticePoint::new(1, 0, 0),
    LatticePoint::new(0, 1, 0),
    LatticePoint::new(0, 0, 1),
];
let bz = brillouin_zone(&lattice);

// Count distinct structures (Burnside's lemma)
let orbits = count_orbits(&PointGroup::c1(), 2); // 2³ = 8 orbits under identity
```

## API

| Type/Function | Description |
|---------------|-------------|
| `LatticePoint` | 3D point in (Z/3Z)³ |
| `SymmetryOp` | Enum: Identity, Inversion, C₃(axis), Mirror(axis) |
| `PointGroup` | Named symmetry group (C₁, C₃ᵢ, Oₕ) |
| `UnitCell` | Basis atoms + lattice vectors |
| `brillouin_zone(&lattice)` | Wigner-Seitz cell in reciprocal space |
| `structure_factor(&cell, &q)` | Diffraction intensity (mod 3) |
| `count_orbits(&group, size)` | Distinct structures via Burnside |

## Architecture Notes

The ternary crystal models the **γ + η = C** conservation principle through group theory:

- **γ (structure)**: the lattice — the periodic arrangement of ternary values
- **η (dynamics)**: symmetry operations that transform the lattice while preserving structure
- **C (conservation)**: the orbit count — the number of genuinely distinct configurations is invariant under relabeling, and Burnside's lemma computes this exactly

The GF(3) arithmetic ensures that all operations stay within the ternary domain — there is no leakage to other number systems. The lattice is closed, complete, and finite per unit cell.

## References

- Weyl, H. (1952). *Symmetry*. Princeton University Press.
- Burns, G. & Glazer, A.M. (2013). *Space Groups for Solid State Scientists* (3rd ed.). Academic Press.
- Lidl, R. & Niederreiter, H. (1997). *Finite Fields* (2nd ed.). Cambridge — GF(3) arithmetic.
- Bravais, A. (1850). *Mémoire sur les systèmes formés par les points distribués régulièrement sur un plan ou dans l'espace*. — Original lattice classification.

## License: MIT

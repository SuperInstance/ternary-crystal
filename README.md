# ternary-crystal

**Crystallography in three-valued space. Point groups, Brillouin zones, and the structure of ternary lattices.**

A crystal is a periodic arrangement of atoms. In real crystallography, atoms sit at real-valued positions and symmetry groups are continuous. But what happens when the lattice is ternary — every coordinate is {-1, 0, +1}?

The answer: the symmetry groups collapse to finite, computable sets. Instead of infinite rotational symmetry, you get C3 (120° rotations). Instead of continuous point groups, you get a handful of discrete operations. The Brillouin zone shrinks to exactly 27 points. Diffraction patterns become exact arithmetic. Everything becomes *computable by hand*.

This crate implements the full crystallographic pipeline for ternary space: lattice points with mod-3 arithmetic, point groups (C1, C3i, Oh), unit cells (simple cubic, body-centered), Brillouin zones, structure factors, and Burnside's lemma for counting distinct crystal structures modulo symmetry.

## What's Inside

- **`LatticePoint`** — a 3D point in Z₃³ with mod-3 dot product and ternary distance
- **`SymmetryOp`** — identity, inversion, C3 rotation, mirror
- **`PointGroup`** — C1, C3i, and Oh (full cubic) symmetry groups
- **`UnitCell`** — simple cubic and body-centered ternary cells
- **`brillouin_zone(lattice)`** — the Wigner-Seitz cell in reciprocal space
- **`structure_factor(cell, q)`** — diffraction intensity at scattering vector q
- **`count_orbits(group, size)`** — Burnside's lemma: how many distinct structures?

## Quick Example

```rust
use ternary_crystal::*;

// Body-centered cubic ternary cell
let cell = UnitCell::body_centered();
let points = cell.generate(1);

// Diffraction pattern
let q = LatticePoint::new(1, 1, 1);
let sf = structure_factor(&cell, &q);
println!("Structure factor at (1,1,1): {}", sf);

// How many distinct structures under Oh symmetry?
let group = PointGroup::ternary_cubic();
let orbits = count_orbits(&group, 3);
println!("Distinct orbits: {}", orbits);
```

## The Deeper Truth

**Ternary crystallography is crystallography stripped to its algebraic bones.** The 230 space groups of real crystals reduce to a small handful of ternary groups. The continuous Brillouin zone becomes 27 discrete points. Structure factors — which in real crystals involve complex exponentials evaluated at real-valued scattering vectors — become simple dot products modulo 3.

Burnside's lemma, which counts the number of distinct colorings of a lattice under symmetry operations, becomes exact and trivial. In real crystallography, this is a hard combinatorial problem. In ternary, you just enumerate: for each symmetry operation, count the fixed points, average. Done.

The structure factor `F(q) = Σ exp(2πi q·r/3)` over atoms at positions r. In ternary, the exponential takes only 3 values (1, ω, ω²) where ω is the cube root of unity. The structure factor is a sum of third roots of unity. This means: F(q) can be 0 (extinction), 1 (constructive), or ω/ω² (partial). Extinctions in real crystals — systematic absences in diffraction patterns — have an exact ternary analog.

**Use cases:**
- **Education** — the simplest possible crystallography, everything computable by hand
- **Ternary lattice design** — enumerate distinct ternary crystal structures
- **Diffraction simulation** — exact diffraction patterns for ternary lattices
- **Symmetry analysis** — point groups and orbit counting for discrete lattices
- **Agent coordination** — crystal symmetry as a metaphor for fleet organization

## See Also

- **ternary-lattice** — lattice structures and tiling
- **ternary-ring** — algebraic ring structure (Z₃ arithmetic)
- **ternary-matrix** — matrix operations over Z₃
- **ternary-sheaf** — sheaf cohomology for data on lattices
- **ternary-topology** — topological analysis of discrete structures
- **ternary-quantum** — quantum mechanics shares crystallographic structure

## Install

```bash
cargo add ternary-crystal
```

## License

MIT

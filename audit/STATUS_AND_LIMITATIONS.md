# Research status and limitations

Snapshot date: 2026-11-01.

The archived calculations are internally consistent for the stated inputs and finite search grids, but the current evidence does not approve a complete classical 128-bit parameter set.

## What has been completed

- The fixed-composition Lee-metric protocol formulation and its implementation checks.
- Exact composition counting and the main correctness and zero-knowledge evidence used by the manuscript.
- The finite attack ledgers recorded in `security/` and `experiments/`.
- The six final `k=60` BKZ calibration runs and the 30-sample PIT summary.
- KMP-style and structured Wagner-style PKP estimates, reported in their native accounting units.

## What the evidence means

The results support the statement that no strategy below the target budget was constructed within the explicitly enumerated models, grids, memory convention, and accounting units. They do not prove a lower bound against every attack. BKZ extrapolation, PKP list-operation estimates, and the finite attack searches must be read with the limitations in the corresponding reports.

## Known gaps

- Some dual, \(\ell_1\)-BDD, algebraic, composition-aware, and ISD-plus-lattice variants are not fully priced.
- Three-layer and higher BJMM-style constructions are not independently priced.
- The PKP multi-list sensitivity values are proxies unless a complete disjoint-block construction and collision analysis is supplied.
- The cost units across all attack families have not been reduced to one uniformly validated classical gate model.
- Full-size BKZ costs remain model extrapolations from finite-dimensional experiments.

For these reasons, use the qualified wording “assessment toward a classical 128-bit target under the stated models” unless the remaining gaps are closed and independently reviewed.

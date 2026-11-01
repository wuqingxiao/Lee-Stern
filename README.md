# Lee-Stern Lee-Metric Identification Protocol

This repository archives the current reproducible materials for the Lee-metric Stern-type identification protocol study.

## Scope

The project studies a corrected fixed-composition formulation of a Stern-type protocol over a Lee metric, together with correctness proofs, implementation checks, and a classical attack assessment using ISD/LMMT, lattice/BKZ, PKP, and related models.

The current evidence supports a **classical 128-bit target assessment under explicitly stated attack models and accounting assumptions**. It does not constitute an unconditional proof that the parameter set has 128-bit security against every possible attack. In particular, some attack families and a common cost unit remain outside the completed assessment, and the full-size BKZ estimates involve model extrapolation.

## Contents

- `paper/`: current manuscript, supplement, preview PDF, and paper-facing evidence notes.
- `protocol/`: reference Python implementation and state-rule tests.
- `counting/`: exact Lee-composition counting checks.
- `security/`: parameter candidates, attack ledgers, verification scripts, and scope statements.
- `experiments/`: reproducible P1/P2 experiment code and recorded outputs.
- `audit/`: independent audit notes and data-consistency checks.

## Reproduction

1. Read `paper/阅读说明.md` and `paper/快速复现.md`.
2. Create a Python 3.11+ environment.
3. Run the scripts in `protocol/`, `counting/`, and `security/` from their respective directories.
4. Run the small-scale BKZ and PKP scripts in `experiments/` only after reviewing their resource requirements.

The recorded results are provided as research evidence and should be independently regenerated before being used as a security claim.

## Status

The repository is an auditable research snapshot dated 2026-11-01. The manuscript deliberately uses qualified language for the 128-bit result and records known limitations in `security/THEORY_AND_SCOPE.md` and `audit/MODEL_AUDIT.md`.

# Changelog

## 0.1.0 - Initial public release candidate

- Riemannian Hager-Zhang line search with standard or approximate Wolfe acceptance.
- Exact retraction derivatives for the sphere, Stiefel, fixed-rank, generalized
  Grassmann and oblique manifolds.
- Nine benchmark figures comparing Manopt's default Armijo line search and RHZLS,
  with normalized convergence quantities and per-iteration cost evaluations.
- Reproducible seeds, disabled gradient stopping, and per-example iteration budgets.
- A single `run_all_experiments` entry point and 30 regression checks.
- Validated with MATLAB R2025a and an external official Manopt 8.0 installation.

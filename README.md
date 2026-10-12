# RHZLS — Riemannian Hager–Zhang Line Search

A MATLAB implementation of the Hager–Zhang line search for Riemannian
optimization, with exact retraction derivatives and examples comparing
RHZLS with Manopt's default Armijo line search.

**Author:** Marco Sutti.

## Installation

Download or clone RHZLS and separately download the official
[Manopt 8.0 release](https://www.manopt.org/downloads/Manopt_8.0.zip).
Extract both packages, then run the following in MATLAB, replacing the
paths with your own:

```matlab
cd('/path/to/RHZLS')
clear functions
startup_rhzls('/path/to/Manopt_8.0')
```

Tested with MATLAB R2025a and Manopt 8.0. The Manopt path can point to the
extracted release root or its inner `manopt` folder. The selection is
retained for the current MATLAB session; `clear functions` resets it.
Manopt is not included in the Git repository.

## Quick start

```matlab
[x, f, info, details] = run_rhzls('sphere');   % Rayleigh quotient
run_rhzls('stiefel');                        % Brockett function
```

Each example compares RHZLS and Manopt Armijo using the same Riemannian
steepest-descent solver, problem, initial point, and retraction.
The four outputs describe RHZLS; both methods' results are in
`details.comparison`. Use `run_rhzls('all')` to run both examples.

## Reproduce the experiments

After setup, one command runs all 58 optimizations, saves three MAT files,
and generates nine PDF figures in `docs/`:

```matlab
results = run_all_experiments;
```

Use `run_all_experiments(false)` to save the results and PDFs without
opening figure windows. To redraw the PDFs from existing MAT results:

```matlab
run('docs/plot_validation.m')
```

Each PDF shows normalized cost errors, gradient norms, and solution errors above, and
cost evaluations per line search, `nf(k)`, below. The benchmarks use fixed
random seeds and disable gradient tolerance stopping to show the plateaus.

| Example | Manifold |
|---|---|
| Rayleigh quotient | Unit sphere |
| Brockett function | Stiefel |
| Positive Rayleigh quotient | Unit sphere |
| Positive Brockett function | Stiefel |
| Orthogonal Procrustes | SO(4) |
| Best rank-2 approximation | Fixed-rank matrices |
| Generalized eigenvalue problem | Generalized Grassmann |
| Thomson tetrahedron | Product of spheres |
| Max-Cut SDP on the cycle C5 | Product of unit circles |

Problem definitions and recorded results are in
[Accuracy examples](docs/ACCURACY_EXAMPLES.md) and
[Additional Manopt examples](docs/MANOPT_CANDIDATES.md).

## Use RHZLS with Manopt

For your own problem, provide Manopt geometry in `problem.M`,
`problem.cost`, `problem.egrad`, and a curve callback
`[y, velocity] = problem.curve(x, d, t)` returning the retracted point and
the **exact derivative of the retraction curve**. Then use:

```matlab
problem.M.retr = problem.curve;
options.linesearch = @rhzls.manopt_linesearch;
options.rhzls = struct('acceptance', 'approximate');
[x, f, info] = steepestdescent(problem, x0, options);
```

RHZLS defaults to `'either'`, accepting standard or approximate Wolfe
conditions; `'approximate'` requires approximate Wolfe acceptance.
A generic vector transport cannot replace the retraction derivative.
See [Advanced usage](docs/USAGE.md) for a complete example and options,
and [Retractions](+rhzls/+retractions/README.md) for the implementations.

## Tests

```matlab
report = test_rhzls;
assert(report.failed == 0);
```

The suite contains 30 checks. See [Validation](docs/VALIDATION.md),
[Manopt baseline](docs/MANOPT_BASELINE.md), and
[Contributing](CONTRIBUTING.md) for further details.

## How to cite

If you use RHZLS for your research, please kindly cite the paper
*Riemannian multigrid line search for low-rank problems*, M. Sutti,
B. Vandereycken, SIAM J. Sci. Comput., 43(3), A1803–A1831, 2021.
[DOI: 10.1137/20M1337430](https://doi.org/10.1137/20M1337430).

```bibtex
@article{SuttiVandereycken2021,
  author  = {Marco Sutti and Bart Vandereycken},
  title   = {Riemannian multigrid line search for low-rank problems},
  journal = {SIAM Journal on Scientific Computing},
  volume  = {43},
  number  = {3},
  pages   = {A1803--A1831},
  year    = {2021},
  doi     = {10.1137/20M1337430}
}
```

Please also identify the RHZLS version or commit used.
[CITATION.cff](CITATION.cff) provides citation metadata.

## License

Copyright (C) 2026 Marco Sutti. Distributed under
[GPL-3.0-or-later](LICENSE). Scientific and third-party credits are in
[NOTICE.md](NOTICE.md).

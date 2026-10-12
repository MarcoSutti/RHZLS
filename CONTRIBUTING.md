# Contributing to RHZLS

RHZLS is maintained by Marco Sutti. Bug reports, reproducible numerical examples,
documentation improvements, and focused pull requests are welcome.

## Report a problem

Include the MATLAB release, operating system, Manopt version, exact command,
random seed, and complete error message. For numerical issues, include a small
problem that reproduces the behavior and state which convergence quantity is
affected. Please remove private data from examples and logs.

## Develop and test a change

Use the official [Manopt 8.0 release](https://www.manopt.org/downloads/Manopt_8.0.zip)
as an external dependency. Do not modify or commit its source files.

Create a branch for one change, then run the regression suite from the RHZLS
directory in MATLAB:

```matlab
clear functions
startup_rhzls('/path/to/manopt')
report = test_rhzls;
assert(report.failed == 0 && report.passed == 30)
```

MATLAB R2025a is the currently validated release. For changes to optimization,
retractions, benchmark definitions, counters, or plotting, also reproduce the
complete experiment suite:

```matlab
results = run_all_experiments(false);
assert(results.totals.runs == 58)
assert(numel(results.paths.pdfs) == 9)
```

The second command saves three MAT files and nine PDFs in `docs/`. Inspect the
affected figures. Keep generated outputs and local dependencies out of pull
requests unless the change explicitly updates published reference results.

## Preserve the scientific comparison

- Compare the methods on the same cost, initial point, manifold, retraction,
  iteration budget, and solver options. Change only the line search.
- Use the official line search selected by Manopt's `steepestdescent` as the
  baseline. Keep Manopt sources unchanged.
- Differentiate the actual retraction curve used to evaluate the cost. A vector
  transport is not a replacement for its derivative.
- Verify new retraction velocities with independent finite differences and the
  tangent constraint. Explain any numerically stable, equivalent slope pairing.
- Keep gradient stopping disabled in the plateau experiments. Show recorded
  histories without artificial padding or clamping positive residuals to `eps`.
- Count actual function and gradient callbacks. Document cache reuse and any
  changes to the definition of a reported quantity.
- Preserve the normalized convergence quantities, transparent evaluation-count
  fills, and legends outside the axes in the published benchmark figures.

## Submit a pull request

Describe the problem, resulting behavior, and validation. State which tests and
experiments you ran, including any failures. Add meaningful regression tests for
algorithmic fixes; a documentation-only change does not need new tests.

The GitHub workflow uses fixed action commits, MATLAB R2025a, and a SHA-256
verified Manopt 8.0 download. It runs the 30 custom RHZLS tests on pushes and pull
requests. Its manual **Run workflow** form also offers the complete experiment
suite. A successful local test run does not establish that the hosted workflow
has run successfully; check its result on GitHub after publication.

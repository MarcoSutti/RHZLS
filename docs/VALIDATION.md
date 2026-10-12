# Current validation

Since 2026-10-12, the package has used the official Manopt snapshot
`6c369671d8834aa28f7e222c7d89c5f216e6953d`, dated 2026-03-03,
whose `manopt_version` reports 8.0.0. All 584 upstream files are included
without changes; provenance and hashes are recorded in `PROVENANCE.txt`
and `UPSTREAM_MANIFEST.json`. The Armijo comparison leaves
`options.linesearch` unset: `steepestdescent` selects its official
default. See [MANOPT_BASELINE.md](MANOPT_BASELINE.md).

The MATLAB R2025a check on 2026-10-12 passed all **30 tests**, with no
failures. The local log is `matlab_official_manopt_tests.log`.
The checks cover the Armijo callback selected by the solver, its original
counters, and rejection of an automatic ng=0 assignment for callbacks
other than the official one.

The full run following the migration to official Manopt, before reducing
the Thomson iteration budget, completed all **58 optimizations**, saved
the three MAT files, and regenerated the nine PDFs. The local log is
`matlab_official_manopt_experiments.log`, which ends with
`OFFICIAL_MANOPT_ALL_EXPERIMENTS_PASS`. The checks confirm that every
history contains all iterations up to its configured maxiter and that
every baseline uses official Armijo, with ng=0 inside the line search.
The totals, including the RHZLS `either` controls that are not plotted,
are nf=300326 and ng=32008. These counters exclude cost and gradient
evaluations performed by the outer solver.

The earlier check of the three candidate examples completed all 18 runs
with official Manopt (`matlab_official_candidate_slopes.log`): 1000
iterations for generalized Grassmann and Thomson, and 60 for Max-Cut.
The Thomson budget was subsequently reduced to 50 for both methods.
The new check of the 18 runs, with Grassmann at 1000 and Max-Cut at 60,
is documented in `matlab_thomson_50.log`; the current MAT files and
Thomson PDF use the new budget.
The same log ends with `THOMSON_50_PASS` and also checks all 58 saved
runs, including those not rerun: each history reaches its configured
maxiter, and each baseline retains the official Armijo callback.
This change did not require rerunning the other suites; the earlier
candidate run with a Thomson budget of 100 is retained in the historical
log `matlab_thomson_100.log`.
The checks verify cost derivatives along retractions using finite
differences and the equivalence of metric and ambient pairings.
For Thomson and Max-Cut, they also verify tangent velocities and unit
column norms. Metric pairings reuse the same egrad and exact velocity,
without additional gradient evaluations.

The `test_rhzls` suite covers line search, retraction derivatives, the
fixed-rank domain, the solver, and the comparison with Manopt Armijo.
The three tests added for evaluation counting check actual cost and
gradient calls, the already available origin, and invalid retraction
attempts.
The `matlab_plateau_*` logs and the checks below document the previous
dependency; they do not validate the new official snapshot.
The check on 2026-10-09 passed all 29 tests and completed the 58 runs
with nine PDFs and exactly nine final windows in interactive mode.
The counter check confirmed that iterations, costs, gradient norms,
steps, and solution distances were identical to the preceding 58
histories. The local log is `matlab_evaluation_histories_validation.log`.
The budgets were subsequently updated by family: 200 for positive
Rayleigh, 850 for positive Brockett, 60 for Procrustes and Max-Cut, and
250 for rank-2 approximation. The base Rayleigh example remains at 250.
The new full run is verified in `matlab_plateau_experiments.log`.
The check after the folder cleanup on 2026-10-08 confirmed 26 passing
tests, no failures, and regeneration of all nine PDFs from the current
MAT files. The local log is `matlab_cleanup_validation.log`; `archive`
remains excluded from the MATLAB path.
The matrix completion example, with and without a cost offset, was
removed. Fixed-rank matrix geometry remains available for rank-2
approximation in the accuracy suite.

The default spectral examples are Brockett/Stiefel and Rayleigh/sphere.
The `run_accuracy_examples` suite adds four families, comparing Armijo,
mixed RHZLS, and RHZLS with approximate Wolfe conditions for three fixed
seeds. Results and criteria are in
[ACCURACY_EXAMPLES.md](ACCURACY_EXAMPLES.md).

The nine PDFs have two panels. The upper panel shows three normalized
quantities: relative cost error, gradient norm relative to its initial
value, and relative solution distance.
The lower panel shows `nf(k)`, the number of cost calls during each line
search. `ng(k)` is saved separately in the MAT files. The counters
exclude outer solver evaluations and diagnostic metrics.
The figures show all recorded iterations, including plateaus.
Gradient tolerance stopping is disabled (`tolgradnorm=0`) in all
experiments. The original data remain in the MAT files.

Three further examples adapted from official Manopt code are described
in [MANOPT_CANDIDATES.md](MANOPT_CANDIDATES.md).

To regenerate the available PDFs without repeating the optimization:

```matlab
startup_rhzls
set(groot, 'defaultFigureVisible', 'off')
run('docs/plot_validation.m')
```

To repeat all 58 runs and generate the nine PDFs with one command:

```matlab
startup_rhzls
results = run_all_experiments;
```

`run_all_experiments(false)` exports without opening windows. The runner
restores the current folder and figure visibility. The local log for the
full migration check is `matlab_official_manopt_experiments.log`; the
subsequent candidate update is in `matlab_thomson_50.log`.

To use an external copy of the official Manopt 8.0 release:

```matlab
cd('C:/path/to/RHZLS')
clear functions
startup_rhzls('C:/path/to/Manopt_8.0')
run_all_experiments
```

The paths refer to the local folders containing the two extracted
packages. Startup also accepts the inner `manopt` subfolder and retains
the selection in subsequent calls without arguments in the same session.
Run `clear functions` before selecting Manopt, because it resets that
selection. The optional check
`report = test_rhzls; assert(report.failed == 0);` can be run after
startup and before the full runner.
The runner recreates problems and results without using callbacks from
pre-existing MAT files that contain paths from the original machine.

The separate check of the numbered Manopt 8.0 release on 2026-10-12
passed all 30 tests and completed all 58 optimizations with nine PDFs
in MATLAB R2025a. It ran in a relocated RHZLS copy, entirely excluding
`external` and selecting the release downloaded from
[Manopt_8.0.zip](https://www.manopt.org/downloads/Manopt_8.0.zip).
The local log `matlab_manopt_8_external_validation.log` ends with
`MANOPT_8_EXTERNAL_PORTABILITY_PASS`: every history contains all
iterations up to maxiter, and every baseline uses the official
`@linesearch` callback, with ng=0. The external check totals are
nf=223351 and ng=26283. The current data distributed with the GitHub
snapshot were not replaced by the external release check.

The earlier `validate_plateau_iterations` check verified all 58 runs:
no iteration was missing up to maxiter. The local log is
`matlab_plateau_iterations.log`. In that check, the nine PDFs had a
764 x 477 pt page size, transparent fills with alpha=0.4, and legends
entirely within the page.

The current check of the nine PDFs confirms a common page size of
1086 x 612 pt, transparency with alpha=0.4, and text entirely within
the page. Visual inspection confirms that the legends are outside the
panels and legible. Selection by root folder, relative path, and
`manopt` subfolder is verified in `matlab_official_manopt_paths.log`;
the three MAT files record the version and commit of the snapshot used.

## First public source check

On 2026-10-12, an archive of the first local Git commit was tested with an
external official Manopt 8.0 installation and no pre-existing result MAT files.
All 30 checks passed, all 58 solver runs reached their configured budgets, and
the nine PDFs were generated. The local log ends with `GITHUB_CLEAN_SOURCE_PASS`.
The workflow YAML was parsed and its action versions and Manopt checksum were
verified. Hosted GitHub Actions has not run yet. `CITATION.cff` passes the
official CFF 1.2 schema, including the preferred 2021 paper citation.

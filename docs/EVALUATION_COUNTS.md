# Function and gradient evaluations per line search

Each comparison PDF has two panels: normalized convergence above, and
the number of objective evaluations per steepest-descent step below.
Both legends are placed outside the data panels, on the right.
The lower panel follows the `nf(k)` comparison used in the earlier
Euclidean experiments, now comparing Manopt Armijo and RHZLS.

## Definitions

- `info(k+1).nf` counts cost callback evaluations during the line search
  producing iterate `k`.
- `info(k+1).ng` counts Euclidean-gradient callback evaluations during that
  same search.
- `info(1)` represents the starting point, iteration zero; both counts are
  zero there.

The cost statistic is `costevals` in both methods; RHZLS also records
`gradevals`. Official Manopt Armijo does not return `gradevals`.
Project postprocessing records `ng = 0` only after validating that the
selected search is the upstream cost-only `@linesearch`. Its raw
statistics and upstream source remain unchanged.
These are actual callback calls, rather than estimates from the number
of inner iterations. Supplying the current cost and directional derivative
avoids evaluating the origin again. Trials rejected by the retraction
domain consume an oracle attempt without calling the cost or gradient.

The counts cover line search only. Initial solver evaluations, gradients
computed by the outer solver, and reference-error diagnostics are excluded.
The supplied benchmarks have separate `cost` and `egrad` callbacks:
Armijo calls only the former during backtracking, while RHZLS calls both
to evaluate the cost and its derivative along the retraction curve.
The counters do not measure runtime or the cost of geometric operations.

The optimizer's caching behavior is retained. Manopt can reuse Armijo's
accepted cost; the current RHZLS adapter leaves accepted-point evaluations
to the outer solver. Those outer evaluations are outside `nf(k)` and `ng(k)`.

## Plotting and reproducibility

`run_all_experiments` runs 58 solver runs, saves all three result MAT files,
and produces nine PDFs. Set its argument to `false` for hidden figures.
Each PDF shows two methods; the accuracy suite additionally saves its
`acceptance='either'` control runs without adding them to the plots.

The lower panel shows `nf(k)` as transparent filled steps over `(k-1,k]`, using Olive
for Armijo and Candy for RHZLS. It excludes iteration zero, ends at each
method's actual final iteration, and shares the convergence panel's
display window. Integer counts use a linear axis. The raw histories remain
complete. Current figures display all recorded iterations.

`ng(k)` remains available in the saved histories. For example:

```matlab
startup_rhzls
data = load('docs/accuracy_examples.mat', 'results');
entry = data.results(1);
armijo = entry.manopt.info;
rhzls = entry.rhzls_approx.info;
sum([armijo.nf])
sum([rhzls.nf])
sum([armijo.ng])
sum([rhzls.ng])
```

Historical MAT files without exact counters require a fresh execution:
`run_all_experiments(false)`. Their Hager–Zhang `evaluations` field counts
oracle attempts and is not substituted for actual function calls.

## Implementation and checks

- `+rhzls/linesearch.m` records actual RHZLS cost and gradient callback calls.
- Official Manopt `linesearch.m` supplies its own `costevals` counter;
  the baseline leaves `options.linesearch` unset.
- `+rhzls/add_evaluation_counts.m` attaches the per-iteration `nf/ng` fields
  and validates the selected line search before assigning Armijo's zero
  line-search gradient count.
- `+rhzls/plot_evaluations.m` draws the lower panel.
- `run_all_experiments.m` saves results and exports each final figure once.

Focused regression tests independently
count callback calls, verify cached-origin exclusion, distinguish invalid
retractions from objective evaluations, and check Armijo backtracking and
the public per-iteration fields. Current checks of the official dependency
are recorded in [VALIDATION.md](VALIDATION.md). The earlier
`matlab_evaluation_counts_tests.log` and `matlab_all_experiments.log`
refer to the previous dependency.
Before the per-example iteration limits were changed on 2026-10-09,
the 58 saved optimization histories were compared with the previous
results: iteration numbers, costs, gradient norms, step sizes, and solution
errors are unchanged. See `matlab_evaluation_histories_validation.log`.

The full migration rerun with official Manopt, before the Thomson limit
was reduced from 1000 iterations, is recorded in
`matlab_official_manopt_experiments.log`: all 58 histories reach their
configured iteration limits, with aggregate line-search counts
`nf = 300326` and `ng = 32008` including the saved accuracy control runs.
The 18 candidate runs were subsequently repeated with Thomson limited to
50 iterations for both methods; see `matlab_thomson_50.log` for the
current candidate histories. That run also validates all 58 saved histories,
without rerunning the other suites. The previous 100-iteration candidate
run is recorded in the historical `matlab_thomson_100.log`.
The preceding aggregate counts belong to
the migration run and do not describe the current shorter suite.
The previous plateau validation is in `matlab_plateau_experiments.log`.
All benchmarks disable gradient stopping and display their complete histories.
The generalized-Grassmann and product-of-spheres slopes project the
already evaluated gradient before pairing it with the same exact
retraction velocity. These equivalent metric pairings reduce cancellation
near stationary points without introducing an additional gradient callback.

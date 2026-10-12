# RHZLS — Riemannian Hager–Zhang Line Search

A MATLAB implementation of the Hager–Zhang line search for Riemannian
optimization, with exact derivatives of retraction curves and reproducible
comparisons against Manopt's default Armijo backtracking.

The package applies the same scalar line-search algorithm to several
manifolds. It includes spectral optimization, orthogonal Procrustes,
fixed-rank approximation, generalized eigenvalue problems, and optimization
over products of spheres. Both methods use the same Manopt
`steepestdescent` solver, problem data, starting point, retraction, and
stopping criteria.

**Author:** Marco Sutti.

**Start here:** install the pinned Manopt dependency, run the quick
examples, or reproduce all nine experiments with `run_all_experiments`.

## Contents

- [Setup and quick start](#setup-and-quick-start)
- [Reproduce the experiments](#reproduce-the-experiments)
- [Examples and saved results](#examples-and-saved-results)
- [Read and customize the plots](#read-and-customize-the-plots)
- [Configure RHZLS](#configure-rhzls)
- [Use RHZLS with your own problem](#use-rhzls-with-your-own-problem)
- [Retractions and their derivatives](#retractions-and-their-derivatives)
- [Tests and package layout](#tests-and-package-layout)
- [How to cite](#how-to-cite)
- [Background and third-party code](#background-and-third-party-code)

## Setup and quick start

The Git repository distributes RHZLS source code and nine reference PDF
figures. **Manopt and generated MAT results are external to the source
checkout.** The supported, pinned dependency is the official
[Manopt 8.0 release ZIP](https://www.manopt.org/downloads/Manopt_8.0.zip),
not a moving checkout of Manopt's `master` branch. Download and extract it
separately. See [MANOPT_BASELINE.md](docs/MANOPT_BASELINE.md) for its
SHA-256 checksum, the comparison contract, and the reference PDF baseline.

Validation uses MATLAB R2025a. Compatibility with older MATLAB releases
has not been validated.

Clone or extract RHZLS, then change MATLAB's current folder to the folder
containing `startup_rhzls.m`. Replace the two paths below with your own:

```matlab
cd('C:/path/to/RHZLS')
clear functions
startup_rhzls('C:/path/to/Manopt_8.0')
set(groot, 'defaultFigureVisible', 'on')

% Run both examples:
[x, f, info, details] = run_rhzls('all');
```

Each quick example displays a convergence figure comparing RHZLS and
Manopt Armijo. For a single example, use `run_rhzls('stiefel')` for Brockett
or `run_rhzls('sphere')` for Rayleigh. Calling `run_rhzls` without arguments
selects Stiefel.
With `'all'`, the outputs are structures with `stiefel` and `sphere`
fields.

To run the regression tests and every experiment with that same Manopt:

```matlab
cd('C:/path/to/RHZLS')
clear functions
startup_rhzls('C:/path/to/Manopt_8.0')
report = test_rhzls;
assert(report.failed == 0);
run_all_experiments
```

The Manopt path may name the extracted release root or its inner `manopt`
folder. After explicit selection, no-argument calls to `startup_rhzls`
reuse it in the same MATLAB session, including those inside the runners
and tests. Run `clear functions` **before** selecting Manopt: doing so
afterwards clears the selection, so you must select the path again.
Use `startup_rhzls` to set the MATLAB path rather than recursively adding
the entire workspace; archived prototypes may have matching names.
The full runner rebuilds all problems and writes fresh MAT results and
PDFs. It does not require any precomputed results.

The external numbered Manopt 8.0 release was validated in MATLAB R2025a
on 2026-10-12 using a relocated RHZLS copy with its entire bundled
`external` folder excluded. All 30 tests passed; the full runner completed
all 58 optimizations through their configured iteration limits and
generated the nine PDFs. The selected baseline was the external release's
official default Armijo callback. This check is separate from the earlier
snapshot validation. Validation logs are generated development artifacts
and are not included in the Git source checkout.

### Complete distribution ZIP

The complete ZIP prepared in the development workspace includes
`external/manopt` and precomputed MAT results in addition to the source
files. That distribution can use `startup_rhzls` without an argument.
Its bundled dependency is pinned to official GitHub commit
`6c369671d8834aa28f7e222c7d89c5f216e6953d` (2026-03-03), retrieved on
2026-10-12. It reports version 8.0.0 and contains changes after the numbered
8.0 release. This ZIP is distinct from GitHub's automatic source ZIP,
which follows the Git repository contents and excludes Manopt and MAT
results. A release download is available only after a release has been
published; the commands above work from the source checkout independently.

## Reproduce the experiments

All commands in this section assume that MATLAB's current folder is
`RHZLS` and `startup_rhzls('C:/path/to/Manopt_8.0')` has already selected
Manopt in the current MATLAB session.

### Repeat every optimization and regenerate the PDFs

For a fresh source checkout, run this first. One command executes all
58 solver runs, writes the three result MAT files, and generates the nine
PDFs:

```matlab
results = run_all_experiments;
```

The runner displays the nine final figures once, after all optimizations
have finished. To export without opening windows, use
`results = run_all_experiments(false);`. It restores the caller's current
folder and previous figure-visibility setting, including when an error
occurs. The returned structure contains `demo`, `accuracy`, `candidates`,
output `paths`, and aggregate line-search counts in `totals`.

### Regenerate all nine PDFs from existing results

After running the experiments, or when using a complete distribution ZIP
with MAT results, this redraws the figures without repeating optimization:

```matlab
startup_rhzls
set(groot, 'defaultFigureVisible', 'on')
run('docs/plot_validation.m')
```

The figures appear in MATLAB and the PDFs are updated in `docs`.
For PDF export without opening figure windows:

```matlab
startup_rhzls
set(groot, 'defaultFigureVisible', 'off')
run('docs/plot_validation.m')
```

Set `defaultFigureVisible` back to `'on'` when you want to see figures
again. This is a MATLAB-wide setting.

`run_rhzls` still does not save its results automatically when called
separately. The two research runners save their MAT files and export their
PDFs by default; pass `false` to `run_accuracy_examples` or
`explore_manopt_examples` to save only the MAT files.

To redraw only one research suite, use `plot_accuracy_examples` or
`plot_manopt_examples`.

## Examples and saved results

| Runner | Example | Manifold |
|---|---|---|
| `run_rhzls('sphere')` | Rayleigh quotient | Unit sphere |
| `run_rhzls('stiefel')` | Brockett function | Stiefel |
| `run_accuracy_examples` | Positive Rayleigh quotient | Unit sphere |
| `run_accuracy_examples` | Positive Brockett function | Stiefel |
| `run_accuracy_examples` | Orthogonal Procrustes | SO(4) |
| `run_accuracy_examples` | Best rank-2 approximation | Fixed-rank matrices |
| `explore_manopt_examples` | Generalized eigenvalue problem | Generalized Grassmann |
| `explore_manopt_examples` | Thomson tetrahedron | Product of four unit spheres |
| `explore_manopt_examples` | Max-Cut SDP on the cycle C5 | Product of five unit circles |

The quick examples use seed 10. The research suites use seeds 1, 2, and 3,
save all runs, and display seed 1 in their PDFs. The accuracy suite performs
36 runs; the Manopt-derived suite performs 18.

| Generated data in `docs/` | PDF output in `docs/` | Figures |
|---|---|---:|
| `demo_results.mat` | `comparison_sphere.pdf`, `comparison_stiefel.pdf` | 2 |
| `accuracy_examples.mat` | `accuracy_*.pdf` | 4 |
| `manopt_candidate_results.mat` | `candidate_*.pdf` | 3 |

Detailed problem definitions, reference solutions, and recorded results
are available in [ACCURACY_EXAMPLES.md](docs/ACCURACY_EXAMPLES.md) and
[MANOPT_CANDIDATES.md](docs/MANOPT_CANDIDATES.md). These supplementary
experiment notes are currently in Italian.

### What is being compared?

RHZLS is the line search used by the outer optimization algorithm.
All comparison experiments use Manopt's Riemannian steepest descent
for both methods.

The Armijo baseline is the official Manopt `steepestdescent` default:
the comparison leaves `options.linesearch` unset and lets the solver
select its own `@linesearch`. The upstream source files are retained
unchanged. Only RHZLS sets a custom line-search callback. See
[MANOPT_BASELINE.md](docs/MANOPT_BASELINE.md) for the pinned source and
baseline parameters. The older RMGLS dependency and the separate copied
Armijo implementation are archived outside the active package.

RHZLS provides two acceptance settings:

- `'either'` accepts standard Wolfe or approximate Wolfe conditions.
  This is the general default and is used by the quick examples.
- `'approximate'` requires the approximate Wolfe conditions. The four
  accuracy PDFs and the three Manopt-derived PDFs use this setting.

The accuracy runner also saves the `'either'` control runs. Its plots
show only Armijo and approximate-Wolfe RHZLS, with no extra control curve.
The baseline in this package is Armijo; the weak-Wolfe baseline in Fig. 5
of the RMGLS paper is a different comparison.

## Read and customize the plots

Each of the nine PDFs has two panels. The upper panel contains three
normalized convergence indicators for each line search:

| Line style | Quantity |
|---|---|
| Solid | Cost gap `abs(f_k - f_star) / abs(f_star)` |
| Dashed | Riemannian gradient norm `norm(g_k) / norm(g_0)` |
| Dotted | Relative solution error, respecting the problem's symmetries |

If `f_star = 0`, the cost gap is divided by the initial gap
`abs(f_0 - f_star)`. Solution errors account for sign, basis, or rotation
ambiguities through the appropriate vector, subspace, or Gram-matrix
comparison.

Manopt uses **Olive circles**; RHZLS uses **Candy triangles**, with filled
markers and black edges. Horizontal reference lines indicate
`sqrt(eps)` and `eps`. An exactly zero error is displayed at `eps` on
the logarithmic axis. Raw costs and gradient norms remain in the saved
histories.

The lower panel shows the number of cost evaluations per line search,
`n_f(k)`, using transparent filled step plots in the same Olive/Candy colors
(`FaceAlpha = 0.4`) so that overlapping regions remain visible. These
integer counts are shown on a linear axis. Both panels share the same
iteration interval, and each method's curve ends at its actual last step.

All recorded iterations are displayed, including the numerical plateaus.
Gradient tolerance stopping is disabled in every benchmark (`tolgradnorm = 0`).

Figures use a common size of approximately half the screen area.
Both legends sit outside the plotting areas on the right; the two panels
retain the same horizontal alignment.
Only PDF figure files are exported.

| Setting | File |
|---|---|
| Palette, fonts, and LaTeX defaults | `options_plot.m` |
| Marker shapes, marker colors, and line styles | `+rhzls/plot_comparison.m` |
| Per-iteration cost evaluation panel | `+rhzls/plot_evaluations.m` |
| Figure dimensions | `+rhzls/figure_position.m` |
| Displayed iteration interval | `+rhzls/plot_window.m` |
| Vector PDF export | `+rhzls/export_pdf.m` |

After changing the plotting settings, rerun `docs/plot_validation.m`.

### Evaluation counts

Every iteration history includes `nf` and `ng`: actual calls to the cost
and Euclidean-gradient callbacks during the line search that produced that
iterate. Both fields are zero at iteration zero. They exclude outer-solver
evaluations, diagnostic metrics, and values already supplied from the
current iterate. A rejected retraction consumes oracle budget but does not
increment either callback count.

For the supplied problems, Armijo evaluates the cost only during its
line search, so `ng = 0` there. Its outer solver still evaluates gradients.
RHZLS evaluates the cost and gradient when computing the retraction-curve
derivative. The gradient counts are saved separately rather than overlaid
on the `nf` plot.

For example, inspect the saved sphere comparison:

```matlab
data = load('docs/demo_results.mat', 'details');
h = data.details.sphere.comparison.rhzls.info;
table([h.iter]', [h.nf]', [h.ng]', ...
    'VariableNames', {'Iteration', 'FunctionEvaluations', 'GradientEvaluations'})
```

The cost counts come directly from `linesearch.costevals`. RHZLS also
records `linesearch.gradevals`. The official Armijo statistics do not
provide `gradevals`; our history postprocessing sets `ng = 0` only after
validating that the selected line search is the upstream cost-only
`@linesearch`. This does not alter Manopt's source or raw statistics.
Hager–Zhang's `linesearch.evaluations` remains its
oracle-attempt budget, which is a different quantity. Function evaluations
alone do not measure runtime: gradient calls and geometric operations also
contribute to the work. See [EVALUATION_COUNTS.md](docs/EVALUATION_COUNTS.md).

Older result MAT files without these counters must be regenerated with
`run_all_experiments(false)` before producing the two-panel figures.

## Configure RHZLS

For the quick examples, problem data belong in `opts.problem`, solver
settings in `opts.solver`, and RHZLS settings in
`opts.solver.linesearch`:

```matlab
opts = struct();
opts.seed = 10;
opts.plot = false;
opts.problem = struct('n', 60);
opts.solver = struct('maxiter', 250, 'tolgradnorm', 1e-14, 'verbosity', 0);
opts.solver.linesearch = struct('acceptance', 'approximate', ...
    'delta', 0.1, 'sigma', 0.9, 'epsilon', 1e-6, 'maxevals', 100);

[x, f, info, details] = run_rhzls('sphere', opts);

% Inspect the final gradient norms of both methods:
comparison = details.comparison;
comparison.manopt.info(end).gradnorm
comparison.rhzls.info(end).gradnorm
```

The four direct outputs describe RHZLS. Both runs, including their
iteration histories, are available in `details.comparison`.
After `run_rhzls('all')`, use `details.sphere.comparison` or
`details.stiefel.comparison`.

The default gradient tolerance is `0`, disabling gradient tolerance
stopping. A positive `tolgradnorm`, as in the customization example above,
enables this criterion for your own runs. The default iteration limits
are 250 for the sphere quick example and 3000 for Stiefel. Both research
runners use a maximum of 1000 iterations, except Positive Rayleigh on the
sphere (200), Positive Brockett on Stiefel (850), Orthogonal Procrustes on
SO(4) (60), best rank-2 approximation (250), Thomson tetrahedron on a
product of spheres (50), and Max-Cut on C5 (60).
Each limit applies to both line searches. Their gradient tolerance is
`0`; settings are specified inside the runner files.

`help rhzls.make_problem` describes the available problem settings,
including a user-supplied symmetric matrix `A`, starting point `x0`,
and Brockett weights.

| RHZLS parameter | Default | Purpose |
|---|---:|---|
| `acceptance` | `'either'` | Standard or approximate Wolfe; use `'approximate'` to require the latter |
| `delta` | 0.1 | Sufficient-decrease parameter |
| `sigma` | 0.9 | Curvature parameter |
| `epsilon` | `1e-6` | Approximate-Wolfe cost tolerance |
| `function_scale` | 1 | Minimum scale for the cost tolerance |
| `alpha0` | 1 | Initial trial step |
| `alphamax` | `1e6` | Upper step limit, further restricted by the retraction domain |
| `rho` | 5 | Bracket expansion factor |
| `gamma` | 0.66 | Required interval contraction |
| `theta` | 0.5 | U3 bisection fraction |
| `maxevals` | 100 | Maximum scalar-oracle evaluations |

Parameters are validated by `rhzls.options`. In particular,
`0 < delta < 1/2` and `delta <= sigma < 1`.

The approximate-Wolfe cost limit is
`phi(0) + epsilon * max(abs(phi(0)), function_scale)`.
It can permit small cost increases. In the solver adapters, the first
trial uses `alpha0`; subsequent trials start from twice the previously
accepted step, subject to the step limit.

For an RHZLS-only run using the package's standalone steepest-descent
solver, set `opts.compare = false`.

## Use RHZLS with your own problem

After installing Manopt 8.0, the following example minimizes a Rayleigh
quotient using Manopt's solver and the RHZLS adapter:

```matlab
startup_rhzls('C:/path/to/Manopt_8.0')

n = 40;
A = diag(linspace(1, 5, n));
[M, curve, maxstep] = rhzls.manifold('sphere', n);

problem = struct('M', M, 'curve', curve, 'maxstep', maxstep);
problem.M.retr = curve;
problem.cost = @(x) x' * A * x;
problem.egrad = @(x) 2 * A * x;
x0 = M.rand();

options = struct('maxiter', 250, 'tolgradnorm', 1e-14, ...
    'minstepsize', 0, 'verbosity', 0);
options.linesearch = @rhzls.manopt_linesearch;
options.rhzls = struct('acceptance', 'approximate');

[x, f, info] = steepestdescent(problem, x0, options);
```

For direct Manopt use, RHZLS parameters go in `options.rhzls`.
The `opts.solver.linesearch` structure above belongs to the
`run_rhzls` interface.

Alternatively, pass the same problem to the standalone solver:

```matlab
solver = struct('maxiter', 250, 'tolgradnorm', 1e-14, 'verbosity', 0);
solver.linesearch = struct('acceptance', 'approximate');
[x, f, info, details] = rhzls.solve(problem, x0, solver);
```

To call the line search directly along a tangent descent direction:

```matlab
g = problem.M.egrad2rgrad(x0, problem.egrad(x0));
d = problem.M.lincomb(x0, -1, g);
[alpha, y, ls] = rhzls.linesearch(problem, x0, d, ...
    struct('acceptance', 'approximate'));
assert(ls.accepted, ls.status);
```

Here `y = R_x0(alpha*d)`. The Manopt adapter returns the step norm
`alpha * problem.M.norm(x0,d)`, as required by Manopt. A failed direct search returns
`accepted = false`, `alpha = 0`, and the original point; the Manopt
adapter raises `rhzls:LineSearchFailed` if no step is accepted.

## Retractions and their derivatives

For a tangent direction `d` at `x`, RHZLS applies Hager–Zhang to

$$
\phi(t) = f(R_x(td)), \qquad
\phi'(t) =
\left\langle \nabla_{\mathrm{E}} f(R_x(td)),
\frac{d}{dt}R_x(td) \right\rangle_F.
$$

The derivative uses the **exact velocity of the chosen retraction curve**,
paired with the Euclidean gradient in the ambient representation.
For generalized Grassmann and products of spheres, an equivalent pairing with the projected
Riemannian gradient and the manifold metric reduces cancellation near
the optimum. The optional `problem.curve_slope(y, egrad, velocity)`
callback receives the already evaluated gradient and that same exact
velocity; it must return the derivative of the cost along the curve.
A generic vector transport cannot replace this velocity.

A custom problem must provide `problem.cost`, `problem.egrad`,
Manopt geometry in `problem.M`, and a two-output curve
`[y, velocity] = problem.curve(x, d, t)`.
Assign the **same function handle** to `problem.M.retr` and
`problem.curve`. An optional `problem.maxstep(x,d)` returns a positive
domain bound or `Inf`.

All implementations are collected in
[`+rhzls/+retractions/`](+rhzls/+retractions/README.md):

| Manifold | Retraction | Implementation |
|---|---|---|
| Unit sphere | Vector normalization | `sphere_curve.m` |
| Stiefel, including SO(n) examples | Polar retraction | `stiefel_curve.m` |
| Fixed-rank matrices | Absil–Oseledets orthographic retraction | `fixedrank_curve.m` |
| Generalized Grassmann | Polar retraction in transformed coordinates | `generalized_grassmann_curve.m` |
| Product of spheres | Columnwise normalization | `oblique_curve.m` |

`rhzls.manifold` constructs sphere, Stiefel, and fixed-rank geometry.
The other two geometries are assembled in `explore_manopt_examples.m`.
Fixed-rank velocities are represented by factors `G,H` with
`velocity = G*H'`; `rhzls.ambient_inner` evaluates the pairing without
forming that dense matrix. `fixedrank_maxstep.m` bounds the smooth
domain of the orthographic retraction.

The older curve entry points in `+rhzls` are compatibility wrappers.
Keep them: saved MAT files can contain function handles that refer to them.

When reusing a problem loaded from a MAT file for further optimization,
restore `problem.M.retr = problem.curve`. MATLAB can deserialize
equivalent anonymous handles as distinct objects. Plot regeneration
does not require this reassignment.

## Tests and package layout

Run the regression suite from the package root, selecting your separate
Manopt 8.0 installation:

```matlab
startup_rhzls('C:/path/to/Manopt_8.0')
report = test_rhzls;
assert(report.failed == 0);
```

The regression suite contains **30 tests**, covering the scalar line search,
exact retraction velocities, pullback derivatives, fixed-rank domain
handling, solver behavior, the Armijo comparison, and callback counts with
cached origins and invalid retractions. All 30 passed with the official
dependency in MATLAB R2025a on 2026-10-12. Current migration checks and
the development logs are described in [VALIDAZIONE.md](docs/VALIDAZIONE.md).
The official-Manopt migration rerun completed all **58 optimizations** and
regenerated the nine PDFs; its log is
`docs/matlab_official_manopt_experiments.log`. Every saved history reaches
its configured iteration limit, and every Armijo baseline uses the actual
official default callback.
After reducing the Thomson limit to 50 iterations, the 18 candidate runs
were repeated; their current log is `docs/matlab_thomson_50.log`. This run
also validates the iteration limits and official Armijo callback in all
58 saved histories; the other suites were not rerun for this limit change.
The preceding 100-iteration candidate run is recorded in the historical
`docs/matlab_thomson_100.log`.
Earlier `matlab_plateau_*` logs document the previous dependency and
remain historical records.
The separate external numbered-release check is recorded in
`docs/matlab_manopt_8_external_validation.log`, ending with
`MANOPT_8_EXTERNAL_PORTABILITY_PASS`. It passed 30 tests and repeated
all 58 runs and nine PDF exports in a relocated package without bundled
Manopt; it did not replace the supplied snapshot's result MAT files.

```text
RHZLS/
├── README.md
├── startup_rhzls.m
├── run_all_experiments.m        Run all suites and export the nine PDFs
├── run_rhzls.m                  Quick spectral examples
├── run_accuracy_examples.m     Four accuracy benchmarks
├── explore_manopt_examples.m   Three Manopt-derived benchmarks
├── plot_accuracy_examples.m    Redraw the accuracy PDFs
├── plot_manopt_examples.m      Redraw the Manopt-derived PDFs
├── options_plot.m              Shared plotting defaults
├── test_rhzls.m
├── +rhzls/                     Common line search and solver adapters
│   └── +retractions/           Retractions and exact curve derivatives
├── tests/                      Regression checks
└── docs/                       Reference PDFs and experiment notes
    └── plot_validation.m       Redraw PDFs after generating MAT results
```

Running the experiments adds three MAT results to `docs/`; these and
local validation logs are ignored by Git. `external/manopt/` is also
ignored and may exist locally or in a complete distribution ZIP. The Git
source checkout uses the separate Manopt 8.0 installation selected above.

In the development workspace, `archive/` is a sibling of `RHZLS/`.
It holds earlier prototypes and intermediate artifacts, with original
paths recorded in `archive/manifest.json`. It is not needed to run
the package and is excluded from the distribution ZIP.

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

For reproducibility, also identify the RHZLS version or Git commit used.
The package is authored by Marco Sutti. [CITATION.cff](CITATION.cff) provides
machine-readable metadata and identifies the paper above as the preferred
citation for GitHub's “Cite this repository” feature.

Please also acknowledge the [Hager–Zhang scalar line-search algorithm](https://doi.org/10.1145/1132973.1132979)
and the [Manopt dependency](https://www.manopt.org/about.html) where appropriate.

## Background and third-party code

The Riemannian formulation follows the approach in Sutti and Vandereycken,
*Riemannian multigrid line search for low-rank problems*. The scalar
bracketing, secant, and bisection core is based on the Hager–Zhang line search
and *Algorithm 851: CG_DESCENT*. The package extends the retraction
derivatives beyond the original fixed-rank setting.

Manopt provides the manifold factories and the common outer solver.
The official Armijo baseline is used directly from the external Manopt
8.0 dependency. Keep its licence and credit files with your Manopt
installation; it is not modified by RHZLS. The dependency version,
download checksum, and provenance of the reference PDFs are documented
in [MANOPT_BASELINE.md](docs/MANOPT_BASELINE.md).
The complete development distribution additionally records snapshot
provenance and upstream licence notices in its `external/manopt` folder.

## License

Copyright (C) 2026 Marco Sutti. RHZLS is distributed under
[GPL-3.0-or-later](LICENSE). Scientific and third-party credits are in
[NOTICE.md](NOTICE.md); citation metadata is in [CITATION.cff](CITATION.cff).

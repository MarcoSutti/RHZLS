# Plotting reference

Generate the result MAT files with `run_all_experiments` before redrawing figures.
Run the commands below from the RHZLS package root after selecting Manopt.
See the [README](../README.md#reproduce-the-experiments) for setup and reproduction.

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
contribute to the work. See [EVALUATION_COUNTS.md](EVALUATION_COUNTS.md).

Older result MAT files without these counters must be regenerated with
`run_all_experiments(false)` before producing the two-panel figures.

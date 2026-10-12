# Official Manopt baseline

The comparisons use the official Manopt source directly. With no
`options.linesearch` override and no problem-level line-search hint,
`steepestdescent` selects its default `@linesearch`: backtracking with the
Armijo sufficient-decrease condition. This differs from the weak-Wolfe
baseline labelled WW in Fig. 5 of the RMGLS paper.

## Pinned dependency for the source repository

The source repository uses the official numbered
[Manopt 8.0 release ZIP](https://www.manopt.org/downloads/Manopt_8.0.zip),
which reports version **8.0.0**, released 2024-07-05. Its SHA-256 is
`b98aab24d9a98634f732edc9ae1d2beb827c83dfc416b1e0732d756574fa7b6a`.
Use this fixed release to reproduce the supported dependency setup;
a moving `master` checkout may change the baseline.

Manopt is installed separately and is not tracked in the RHZLS Git
repository. The repository also excludes generated MAT histories and
local validation logs. It retains nine reference PDFs, described below.

Extract the Manopt release and select it from the RHZLS package root:

```matlab
cd('C:/path/to/RHZLS')
clear functions
startup_rhzls('C:/path/to/Manopt_8.0')
run_all_experiments
```

Replace the paths with your local folders. The Manopt path can identify
the extracted release root or its inner `manopt` folder. The startup
checks that required factories and solvers resolve within that selection.
Subsequent no-argument startup calls retain the selection in the same
MATLAB session. Clear functions before making the explicit selection;
clearing them afterwards resets it. For an optional preliminary check,
run `report = test_rhzls; assert(report.failed == 0);` after startup.

The full runner recreates the problems and writes new MAT histories and
PDFs for the selected Manopt. A fresh source checkout must run it before
using `docs/plot_validation.m`, because the generated MAT histories are
not committed. No Manopt source files are modified.

### Validation of the numbered 8.0 release

On 2026-10-12 the official release ZIP was downloaded independently from
the link above and tested in MATLAB R2025a.
The check used a relocated copy of RHZLS with its entire `external`
directory excluded, selecting only the independently extracted release.
All 30 regression tests passed; all 58 optimization runs reached their
configured iteration limits; nine PDFs were exported. Every Armijo run
used the selected release's actual default `@linesearch` and recorded
`ng = 0`. The log is `matlab_manopt_8_external_validation.log`, ending with
`MANOPT_8_EXTERNAL_PORTABILITY_PASS`.

The validation log is a local development artifact, excluded from Git.

## Reference PDFs and complete distribution ZIP

The nine reference PDFs committed in the source repository were generated
with the official GitHub snapshot at commit
[`6c369671d8834aa28f7e222c7d89c5f216e6953d`](https://github.com/NicolasBoumal/manopt/tree/6c369671d8834aa28f7e222c7d89c5f216e6953d),
dated 2026-03-03 and retrieved on 2026-10-12. It reports 8.0.0 but includes
source changes after the numbered 8.0 release. It is a fixed snapshot,
not a separately numbered newer release.

The complete ZIP prepared in the development workspace additionally
contains that snapshot in `external/manopt`, the corresponding MAT
histories, and local validation records. These files are absent from
GitHub's automatic source ZIP. A complete release attachment is available
only after one is published; it is not required for the setup above.

The 584 upstream snapshot files were copied unchanged. The complete
distribution adds `PROVENANCE.txt` and `UPSTREAM_MANIFEST.json` alongside
them; the manifest records file hashes for integrity checks. Upstream
licence and credit files remain at the snapshot root: `LICENSE.txt`,
`COPYING.txt`, `CREDITS.txt`, and `CLA.txt`.

Rerunning with Manopt 8.0 generates results for that exact release rather
than reproducing the snapshot's floating-point trajectory byte for byte.
Trajectories can also differ across MATLAB releases and platforms.

## Comparison contract

The Armijo run leaves `options.linesearch` unset. RHZLS sets
`options.linesearch = @rhzls.manopt_linesearch`. Both runs receive the same
problem, initial point, exact retraction implementation, iteration limit,
and solver stopping settings. Gradient stopping is disabled
(`tolgradnorm = 0`) so the plots show the recorded numerical plateaus.
The shared retractions are supplied by the RHZLS package; this does not
modify the upstream manifold factories.

The default parameters used by the Manopt 8.0 release and the reference
snapshot's upstream
[linesearch.m](https://github.com/NicolasBoumal/manopt/blob/6c369671d8834aa28f7e222c7d89c5f216e6953d/manopt/solvers/linesearch/linesearch.m)
are:

| Parameter | Value |
| --- | ---: |
| `ls_contraction_factor` | `0.5` |
| `ls_optimism` | `2` |
| `ls_suff_decr` | `2^-13` |
| `ls_max_steps` | `25` |
| `ls_initial_stepsize` | `1` |

The benchmark runners leave these parameters at their upstream defaults.
The [official solver documentation](https://www.manopt.org/solvers.html)
also identifies backtracking as the default for steepest descent.

## Evaluation counts

`nf(k)` is read from the upstream `info(k+1).linesearch.costevals`, without
wrapping or replacing Armijo. The official statistics have no `gradevals`
field. For this validated `@linesearch`, which calls only the cost during
backtracking, project postprocessing records `ng(k) = 0` in the iteration
history. RHZLS records its own cost and gradient callback counts.
Outer-solver gradient evaluations are excluded from both line-search
counts; see [EVALUATION_COUNTS.md](EVALUATION_COUNTS.md).

## Previous baseline

The RMGLS bundle previously included Manopt 5.0.0 with local changes,
including an extra halving of an accepted Armijo step. The project then
used a separately copied Armijo implementation with that halving removed.
Both the old dependency and copied baseline are preserved in
`archive/history/2026-10-12_before_official_manopt`, outside the active
MATLAB path. The active comparisons now call the official default directly.

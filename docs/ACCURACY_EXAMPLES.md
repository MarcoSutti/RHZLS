# Accuracy examples

The suite comprises four problem families with three fixed initializations
(seeds 1, 2, 3), for a total of 36 runs. The current results are saved in
`accuracy_examples.mat`. Since October 12, 2026, the comparison has used the
official Manopt snapshot described in [MANOPT_BASELINE.md](MANOPT_BASELINE.md),
allowing `steepestdescent` to select its own default Armijo line search
directly. The migration checks are documented in
[VALIDATION.md](VALIDATION.md).
The outer solver is Manopt's `steepestdescent`, with the same problem,
initial point, retraction, and `tolgradnorm=0` for all three methods.
Since October 9, 2026, `maxiter` has been set to 200 for positive Rayleigh,
850 for positive Brockett, 60 for Procrustes, and 250 for rank-2 approximation.

## The distinction between line-search acceptance modes

The earlier version uses `acceptance='either'`: it terminates when either
the ordinary Wolfe conditions or the approximate Wolfe conditions hold.
Near a minimum, the cost can round to the same value at different points,
and the ordinary sufficient-decrease condition can appear satisfied
numerically even for excessively large steps. The new scalar test
reproduces this situation with a quadratic cost whose base value is `1e20`:
the approximate mode rejects the step that overshoots the minimum.

The explicit option `acceptance='approximate'` was added to require both
derivative bounds as well as the cost bound. The mixed mode remains the
default and is not presented as equivalent to the approximate mode.
Hager and Zhang describe modes V1, V2, and V3 (Wolfe, approximate Wolfe,
or a permanent switch from one to the other) in Section 3, page 121, of
[Algorithm 851: CG_DESCENT](https://people.clas.ufl.edu/hager/files/cg_code.pdf).
This project now implements an explicit choice between the existing mixed
mode and acceptance based only on approximate Wolfe conditions; it does
not implement the adaptive V3 policy or the complete CG_DESCENT conjugate
gradient solver.

## Problems without artificial cost shifts

1. **Positive Rayleigh on the sphere S^11**: `f(x)=x'*D*x`, with
   `D=diag(linspace(1,5,12))`, minimum 1, and minimizers +/-e1.
2. **Positive Brockett on Stiefel St(12,2)**: `trace(X'*D*X*N)`,
   with the same D and `N=diag([2,1])`. The first two eigenvectors minimize
   the cost, up to column signs.
3. **Orthogonal Procrustes on SO(4)**: `f(Q)=0.5*norm(Q-A,'fro')^2`,
   with `A=diag(1:4)`. The minimum is 7, attained at Q=I. The square Stiefel
   factory is used, initialized with positive determinant; the polar
   retraction preserves this component of the orthogonal matrices.
4. **Best rank-2 approximation**, on the manifold of 12x12 fixed-rank
   matrices: `f(X)=0.5*norm(X-A,'fro')^2`, with
   `A=diag(linspace(2,0.2,12))`. The solution is the truncated SVD; the
   positive minimum is half the sum of the squared discarded singular values.

The following table reports the final gradient norms from the 36 runs with
the official Manopt snapshot. The values are confirmed by the current
results and the log `matlab_official_manopt_experiments.log`.

| Example | Armijo gradient norm, range over 3 seeds | RHZLS approx. gradient norm, range over 3 seeds | RHZLS approx. iterations |
|---|---:|---:|---|
| Positive Rayleigh | 2.981e-8 -- 1.025e-7 | 5.070e-18 -- 1.987e-15 | 200 for every seed |
| Positive Brockett | 4.216e-8 -- 4.230e-7 | 2.311e-17 -- 1.233e-15 | 850 for every seed |
| Orthogonal Procrustes | 3.026e-8 -- 8.843e-8 | 5.636e-32 -- 5.163e-27 | 60 for every seed |
| Rank-2 approximation | 2.441e-8 -- 3.623e-8 | 4.660e-15 -- 3.823e-13 | 250 for every seed |

All methods complete their specified iteration budgets: gradient tolerance
does not trigger termination. Small norms are computed residuals, not a
measure of the number of accurate digits: in structured problems, some
coordinates can continue to decrease below eps. No artificial plateaus
are imposed, and no iterations are added to the histories.

Procrustes and fixed-rank approximation are recommended for an illustrative
figure: both have a natural positive minimum and a known solution.

## Reproduction

```matlab
cd('C:/path/to/RHZLS')
startup_rhzls
set(groot, 'defaultFigureVisible', 'on')
results = run_accuracy_examples;
```

To regenerate only the new figures from the saved results:

```matlab
plot_accuracy_examples
```

For each family, a single `accuracy_*.pdf` shows seed 1 and three metrics
for Armijo and approximate RHZLS (curves for the `either` mode are excluded
from the plots): relative cost error, gradient norm divided by the common
initial gradient norm, and relative distance to the solution. Separate
`_control` figures are not generated. The other two seeds remain available
in the results and the complete log.
Cost errors are divided by `abs(f_star)`; when the minimum is zero,
`abs(f_0-f_star)` is used. The original MAT data remain absolute.
Values that round to zero are displayed at eps. Distances and gradients
are evaluated independently of cost subtraction. Curves end at the actual
termination of each run.

The displayed range includes all iterations: 200 for the sphere, 850 for
Stiefel, 60 for Procrustes, and 250 for fixed rank.

To select this mode in an existing example:

```matlab
opts.solver.linesearch = struct('acceptance','approximate');
run_rhzls('sphere', opts);
```

The regression suite includes a check for an excessive step that is
acceptable only because of cost rounding. The current test results with
official Manopt are reported in [VALIDATION.md](VALIDATION.md).

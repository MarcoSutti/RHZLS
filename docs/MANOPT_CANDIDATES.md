# Additional examples from Manopt

The script `explore_manopt_examples.m` adapts three problems from the
official [Manopt examples](https://github.com/NicolasBoumal/manopt/tree/master/examples)
to study the comparison between Armijo and RHZLS with approximate Wolfe
conditions. These experiments form a separate suite from the existing
benchmarks. The cost functions have naturally nonzero minima: no cost
offset is added. Differences in accuracy must be measured, without
assuming that they will appear for every problem or starting point.

## Reproducible protocol

Each problem uses seeds `1:3` with `rng(seed, 'twister')`.
For each seed, both methods receive the same problem, starting point,
and retraction. Both call Manopt's `steepestdescent`, with
`maxiter = 1000` for generalized Grassmann, 50 for Thomson, and 60 for
Max-Cut on C5, `tolgradnorm = 0`, `minstepsize = 0`, and the same
remaining options.

The only difference is the selected line search:

- Official Manopt: no override of `options.linesearch`, so
  `steepestdescent` selects its default Armijo line search, `@linesearch`;
- `rhzls.manopt_linesearch`, with `acceptance = 'approximate'`: RHZLS with
  acceptance through approximate Wolfe conditions.

The provenance of the Armijo baseline is described in `MANOPT_BASELINE.md`.
The official examples also call other solvers: here, we adapt their cost
functions and geometries, using steepest descent to isolate the
contribution of the line search. A total of 18 optimizations are run.

The script also checks the retraction velocity and the derivative of the
cost along the curve against centered finite differences (`t = 0.15`,
`h = 1e-6`, required relative error below `1e-7`), and checks that the
velocity at `t = 0` matches the tangent direction.

The plotted quantities are always normalized:

\[
\frac{|f_k-f_\star|}{|f_\star|},\qquad
\frac{\|\operatorname{grad} f(x_k)\|}{\|\operatorname{grad} f(x_0)\|},
\qquad e_{\mathrm{solution}}(x_k).
\]

Gradient norms and absolute cost errors remain available in the raw
results. The solution error is already relative and accounts for the
symmetries of each problem.

## 1. Dominant eigenspace on generalized Grassmann

Source: [generalized_eigenvalue_computation.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/generalized_eigenvalue_computation.m).

We minimize

\[
f(X)=-\operatorname{tr}(X^T A X),\qquad X^T B X=I_2,
\]

on `grassmanngeneralizedfactory(12, 2, B)`. In the controlled instance,
`B = diag(linspace(1,2,12))`, `L = chol(B,'lower')`,
`lambda = linspace(1,5,12)`, and `A = L*diag(lambda)*L'`.
The whitened coordinates are `Q = L'*X`: the optimal subspace is spanned
by the last two standard basis vectors. Therefore,

\[
f_\star=-(\lambda_{11}+\lambda_{12})=-\frac{106}{11}.
\]

The starting point is obtained from the QR factorization of a Gaussian
`12 x 2` matrix. The solution error is

\[
e_{\mathrm{solution}}(X)=\frac{\|QQ^T-P_\star\|_F}{\sqrt{2}},
\]

where `Pstar` is the projector onto the last two coordinates. This avoids
penalizing a rotation of the basis within the same subspace.

The polar retraction is used in whitened coordinates:

\[
R_X(tD)=L^{-T}\operatorname{polar}(L^T(X+tD)),\qquad
\frac{d}{dt}R_X(tD)=L^{-T}\frac{d}{dt}
\operatorname{polar}(L^T(X+tD)).
\]

The polar derivative is computed by `rhzls.retractions.stiefel_curve`.
The script explicitly assigns this same retraction to both methods,
so the comparison uses identical retractions.

## 2. Thomson problem: a tetrahedron on a product of spheres

Source: [thomson_problem.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/thomson_problem.m).

Four points on `S^2` are represented by the columns of
`X` on `obliquefactory(3,4)`. The official normalized cost is retained:

\[
f(X)=\frac{1}{16}\sum_{i<j}\frac{1}{1-x_i^T x_j}.
\]

The global minimum is a regular tetrahedron, with off-diagonal inner
products equal to `-1/3`, and

\[
f_\star=\frac{6}{16}\frac{3}{4}=\frac{9}{32}.
\]

Indeed, `norm(sum(X,2))^2 >= 0` implies that the average pairwise inner
product is at least `-1/3`. The convexity of `1/(1-s)` gives the lower
bound above, attained by the tetrahedron.

The official example initially proposes a finite-difference gradient,
which is unsuitable for this accuracy comparison. Here, the ambient
gradient is analytical:

\[
\nabla f(X)=\frac{XW}{16},\qquad
W_{ij}=(1-x_i^T x_j)^{-2}\ (i\ne j),\qquad W_{ii}=0.
\]

The starting point is random on the manifold. The solution error uses
the reference Gram matrix `Gstar = (4/3)*eye(4)-ones(4)/3`:

\[
e_{\mathrm{solution}}(X)=
\frac{\|X^T X-G_\star\|_F}{\|G_\star\|_F}.
\]

This error is invariant under rotations and, for the regular tetrahedron,
under vertex permutations. Coincident points give an infinite cost;
the script explicitly handles this case.

## 3. Max-Cut relaxation on the cycle C5

Sources: [maxcut.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/maxcut.m)
and [elliptope_SDP.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/elliptope_SDP.m).

We consider the Laplacian `L` of the unweighted cycle with five vertices,
and five unit vectors in the plane, stored as the columns of
`X` on `obliquefactory(2,5)`. The cost of the relaxation is

\[
f(X)=-\frac14\operatorname{tr}(XLX^T),\qquad
\nabla f(X)=-\frac12 XL.
\]

This is the transpose of the representation with unit rows used in the
official Max-Cut example. Here, the rank remains fixed at two: the rank
increases and rounding of the full program are not applied.

The reference has columns

\[
x_{\star,j}=\begin{bmatrix}
\cos(4\pi j/5)\\ \sin(4\pi j/5)
\end{bmatrix},\qquad j=0,\ldots,4,
\]

and minimum relaxation cost

\[
f_\star=-\frac52(1+\cos(\pi/5))\simeq-4.522542486.
\]

The spectral bound follows from `lambda_max(L) = 2+2*cos(pi/5)`
and `trace(X*X') = 5`; the pentagram configuration attains it.
The starting point is random. The solution error is the relative Gram
matrix difference from `Xstar'*Xstar`, so orthogonal transformations in
the plane do not change the indicator.

## Retraction and velocity for the two oblique problems

For each column, with `z_i = x_i+t*d_i` and `s_i = norm(z_i)`,

\[
y_i=\frac{z_i}{s_i},\qquad
\dot y_i=\frac{d_i-y_i(y_i^T d_i)}{s_i}.
\]

The function `oblique_curve` in `+rhzls/+retractions`, also called by the
local wrapper in `explore_manopt_examples.m`, evaluates these formulas.
Normalization is independent for each column: the entire matrix is not
normalized by a single norm.

## Candidates for future evaluation

- [positive_definite_intrinsic_mean.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/positive_definite_intrinsic_mean.m):
  intrinsic mean of SPD matrices. With commuting, well-conditioned data,
  an exact reference is available through the mean of the logarithms.
  The local retraction is `X+t*D+0.5*t^2*D*(X\D)`, with velocity
  `D+t*D*(X\D)`. The official code warns that gradient evaluation can
  limit accuracy when the noise is large: this limitation should be
  distinguished from that of the line search.
- [doubly_stochastic_denoising.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/doubly_stochastic_denoising.m):
  projection in the Frobenius norm onto positive doubly stochastic
  matrices. This would require the derivative of Sinkhorn scaling.
  The Manopt factory included in the project uses PCG with relative
  tolerance `1e-6` for the gradient and projection, together with
  iterative scaling and final clipping. These numerical errors can
  dominate those of the line search at very high accuracy. This
  problem is not included in the pilot suite.

## Running the experiments and results

From MATLAB:

```matlab
cd('C:/path/to/RHZLS')
startup_rhzls;
results = explore_manopt_examples();
plot_manopt_examples(results);
```

Current results are saved in `manopt_candidate_results.mat`.
Since October 12, 2026, the official snapshot described in
`MANOPT_BASELINE.md` has been used; the migration checks are reported in
`VALIDATION.md`. The complete histories also include any cases that do
not converge within the budget and explicit solver errors.
The validation with official Manopt and the updated Thomson budget,
documented in `matlab_thomson_50.log`, successfully ran all
18 optimizations and derivative checks. Both methods complete
1000 iterations for Grassmann, 50 for Thomson, and 60 for Max-Cut.
The gradient tolerance is disabled. The previous experiment with
100 iterations is documented in the historical log
`matlab_thomson_100.log`.

The following table reports current results with official Manopt.

| Problem | Seed | Final Armijo norm | RHZLS iterations | Final RHZLS norm | RHZLS relative solution error |
|---|---:|---:|---:|---:|---:|
| Generalized Grassmann | 1 | `3.686e-8` | 1000 | `1.374e-92` | `1.342e-92` |
| Generalized Grassmann | 2 | `2.422e-7` | 1000 | `4.779e-15` | `3.945e-16` |
| Generalized Grassmann | 3 | `1.391e-8` | 1000 | `1.633e-15` | `3.729e-16` |
| Thomson, tetrahedron | 1 | `5.408e-9` | 50 | `1.719e-17` | `1.799e-16` |
| Thomson, tetrahedron | 2 | `3.991e-9` | 50 | `1.252e-17` | `1.862e-16` |
| Thomson, tetrahedron | 3 | `2.667e-11` | 50 | `8.845e-18` | `5.888e-17` |
| Max-Cut, cycle C5 | 1 | `3.191e-8` | 60 | `4.839e-16` | `2.243e-16` |
| Max-Cut, cycle C5 | 2 | `2.657e-8` | 60 | `6.131e-16` | `3.621e-16` |
| Max-Cut, cycle C5 | 3 | `1.757e-8` | 60 | `1.205e-15` | `4.844e-16` |

For Grassmann, the derivative along the same retraction is evaluated as
`M.inner(y, M.egrad2rgrad(y, egrad), velocity)`. The exact velocity is
horizontal: this is mathematically equivalent to the ambient pairing,
but reduces cancellation of the normal component near the minimum.
The projection reuses egrad, without a second gradient call.
Finite differences and the comparison between both pairings are checked
for all seeds. Exactly representable zero coordinates can produce
residuals below eps; this does not imply general accuracy beyond machine
precision.

For Thomson and Max-Cut, we also use
`M.inner(y, M.egrad2rgrad(y, egrad), velocity)`. On products of spheres,
the exact velocity is tangent to each column: removing the normal
component of the gradient before pairing preserves the derivative and
reduces cancellation near a stationary point. Egrad is reused, so no
additional gradient calls are made. The script checks the pairing
against the ambient pairing and finite differences of the cost, as
well as velocity tangency and unit column norms.

The final Armijo accuracy is not always `1e-8`: for Thomson, seed 3,
it reaches `2.667e-11`. The observed accuracy difference in these cases
is therefore documented by the results for all three seeds, rather
than assumed to be universal. Near the minimum, the cost difference
can round to zero even when the gradient and solution error remain
distinguishable: the comparison uses all three indicators.

`plot_manopt_examples.m` produces a single figure per problem, using
seed 1 and the saved histories. Each figure contains the six curves of
the three indicators for Armijo and RHZLS, with colors, markers, and
dimensions shared with the rest of the project. The PDF files are
`candidate_generalized_eigen.pdf`, `candidate_thomson_tetrahedron.pdf`,
and `candidate_maxcut_cycle5.pdf`. The figures show all recorded
iterations, without cropping the plateaus.

To regenerate and display only these PDFs, without rerunning the
optimizations:

```matlab
cd('C:/path/to/RHZLS')
startup_rhzls;
set(groot, 'defaultFigureVisible', 'on')
plot_manopt_examples;
```

The retraction and derivative implementations are collected in
`+rhzls/+retractions`: `generalized_grassmann_curve.m` and `oblique_curve.m`
call the polar derivative in `stiefel_curve.m` where appropriate.
To reuse a `problem` loaded from a MAT file in a new optimization,
reassign `problem.M.retr = problem.curve`: MATLAB can reconstruct two
anonymous handles that are numerically equivalent but not identical
when loading.

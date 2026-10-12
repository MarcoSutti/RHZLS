# Advanced usage

Install and select Manopt as described in the [README](../README.md#installation).
Run the commands below from the RHZLS package root.
For figure settings and evaluation counters, see [Plotting reference](PLOTTING.md).

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
[`+rhzls/+retractions/`](../+rhzls/+retractions/README.md):

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

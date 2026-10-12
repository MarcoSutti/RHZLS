# Retractions and exact curve derivatives

This folder contains the retraction implementations and their derivatives
along the curve `t -> R_x(t*d)` used by RHZLS.

| File | Manifold | Second output |
|---|---|---|
| `sphere_curve.m` | Unit sphere | Derivative of vector normalization |
| `stiefel_curve.m` | Stiefel, including SO(n) examples | Polar derivative obtained through a Sylvester equation |
| `fixedrank_curve.m` | Fixed-rank matrices | Factors `G,H` of the ambient velocity `G*H'` |
| `generalized_grassmann_curve.m` | Generalized Grassmann, `X'*B*X = I` | Polar derivative in transformed coordinates; requires `L = chol(B,'lower')` |
| `oblique_curve.m` | Product of spheres | Derivative of columnwise normalization |

Call these functions through their MATLAB namespace. There is no need
to add this folder separately to the path:

```matlab
[y, velocity] = rhzls.retractions.sphere_curve(x, d, t);
[y, velocity] = rhzls.retractions.generalized_grassmann_curve(x, d, t, L);
```

The curve's second output is its exact ambient velocity. A generic vector
transport is not a substitute for this derivative. The same function
handle must be assigned to `problem.M.retr` and `problem.curve`.

The older functions directly in `+rhzls`, and the local curve functions
in `explore_manopt_examples.m`, delegate to these implementations.
They preserve compatibility with existing scripts and function handles
stored in the supplied MAT files.

See the [main README](../../README.md#retractions-and-their-derivatives)
for the algorithmic contract and an example of using RHZLS with a custom
problem.

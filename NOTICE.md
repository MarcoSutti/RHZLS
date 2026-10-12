# Copyright, licence and scientific credits

Copyright (C) 2026 Marco Sutti.
RHZLS is licensed under the GNU General Public License, version 3 or later
(SPDX: GPL-3.0-or-later). See [LICENSE](LICENSE).

Marco Sutti is the author and maintainer of this software package. Authorship
of the package is distinct from the authorship of the referenced algorithms,
research papers and external dependencies.

## Scientific background

The Riemannian formulation is based on Marco Sutti and Bart Vandereycken,
*Riemannian Multigrid Line Search for Low-Rank Problems*, SIAM Journal on
Scientific Computing 43(3), A1803-A1831 (2021),
https://doi.org/10.1137/20M1337430.
The scalar line-search strategy is based on William W. Hager and Hongchao
Zhang, *Algorithm 851: CG_DESCENT, a Conjugate Gradient Method With Guaranteed
Descent*, ACM Transactions on Mathematical Software 32(1) (2006),
https://doi.org/10.1145/1132973.1132979.
RHZLS implements line search and steepest-descent experiments, not the full
CG_DESCENT conjugate-gradient solver or the RMGLS multigrid algorithm.

## External software and adapted examples

Manopt is an external dependency, copyright Nicolas Boumal and distributed
under GPL version 3 or later, with additional notices for its third-party
components. Retain the upstream notices when redistributing Manopt.
Upstream project: https://www.manopt.org/ and
https://github.com/NicolasBoumal/manopt.

The generalized-eigenvalue, Thomson and Max-Cut examples in
`explore_manopt_examples.m` are adapted from Manopt examples. Their upstream
credits include Bamdev Mishra (generalized eigenvalues) and Nicolas Boumal
(Thomson and Max-Cut). The local instances, gradients, reference metrics and
matched retraction curves are specified in this package. See
[MANOPT_BASELINE.md](docs/MANOPT_BASELINE.md) for the pinned dependency.

The source repository contains RHZLS and its reference figures. The Manopt
source tree is installed separately and is excluded from Git tracking.

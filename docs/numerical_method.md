# Numerical Method

## Spatial and temporal discretization
The code discretizes space on a uniform 1D grid and advances time in fixed steps.

## Crank–Nicolson scheme
Time evolution uses the Crank–Nicolson method, forming linear systems of the form

\[
A\,\psi^{n+1} = B\,\psi^n.
\]

The implementation builds these matrices from the Laplacian stencil and potential contribution.

## Boundary conditions
- **Dirichlet:** interior points are solved with endpoint values fixed.
- **Periodic:** coupling terms connect first and last grid points.

## Linear solvers
For tridiagonal Dirichlet systems, the script includes a Thomas algorithm implementation and also supports dense solving for comparison. The periodic case is solved with dense linear algebra in the current script.

## Diagnostics present in the repository
The organized outputs include diagnostics already produced by the original work:
- Norm conservation checks,
- Convergence/error comparisons,
- Computational cost comparisons,
- Reflection/transmission plots.

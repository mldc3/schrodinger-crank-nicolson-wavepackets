# Numerical Method (Implementation-Focused)

This document explains how the implementation in `src/schrodinger_crank_nicolson.py` is organized to evolve 1D quantum wave packets with Crank–Nicolson, compare linear solvers, and produce diagnostics/figures.

## 1) Spatial grid construction

The code defines a 1D uniform grid with `numpy.linspace`:

$$
x_j = x_{\min} + j\,\Delta x, \quad j=0,\dots,N-1.
$$

From this array, the spacing is obtained numerically (`dx = x[1] - x[0]`). Most operators are built from this constant spacing.

## 2) Time grid

Time is advanced with a fixed step `dt` for a configured number of iterations `steps`:

$$
t^n = n\,\Delta t, \quad n=0,1,2,\dots
$$

The implementation loops over time levels and updates the complex state at each iteration.

## 3) Complex wave-function array

The wave function is represented as a complex NumPy array (`dtype=complex`), storing

$$
\psi_j^n \approx \psi(x_j,t^n).
$$

Using a complex array is essential because both phase and amplitude are needed for Schrödinger dynamics.

## 4) Potential array

Each physical scenario defines a real-valued `V(x)` sampled on the same spatial grid:

$$
V_j = V(x_j).
$$

The script includes free propagation, confined cases, harmonic confinement, and localized scattering potentials (finite well / rectangular barrier).

## 5) Gaussian packet initialization

Initial conditions are built from a Gaussian envelope times a plane-wave phase:

$$
\psi(x,0) = \exp\!\left[-\frac{(x-x_0)^2}{2\sigma^2}\right] \exp(i k_0 x).
$$

Implementation parameters control center (`x0`), width (`sigma`), and carrier wave number (`k0`).

## 6) Normalization

After initialization, the state is normalized numerically so that probability starts at one:

$$
\psi \leftarrow \frac{\psi}{\sqrt{\int |\psi|^2 dx}}.
$$

In code, this is evaluated with a grid-based integral (`np.trapz`/`np.trapezoid` style usage).

## 7) Finite-difference Laplacian

The second derivative in the Hamiltonian is approximated with a centered stencil:

$$
\frac{\partial^2 \psi}{\partial x^2}\Big|_{x_j}
\approx
\frac{\psi_{j+1}-2\psi_j+\psi_{j-1}}{\Delta x^2}.
$$

This produces nearest-neighbor couplings and therefore a tridiagonal structure for Dirichlet interior systems.

## 8) Hamiltonian matrix assembly

The discrete Hamiltonian has the form

$$
H = -D_{xx} + \operatorname{diag}(V),
$$

with entries (Dirichlet interior form)

$$
H_{j,j}=\frac{2}{\Delta x^2}+V_j, \qquad
H_{j,j\pm1}=-\frac{1}{\Delta x^2}.
$$

The implementation constructs these diagonals explicitly (dense or sparse-style paths depending on section).

## 9) Crank–Nicolson matrices A and B

Using

$$
i\,\partial_t\psi = H\psi,
$$

Crank–Nicolson gives

$$
\left(I+\frac{i\Delta t}{2}H\right)\psi^{n+1}
=
\left(I-\frac{i\Delta t}{2}H\right)\psi^n.
$$

Define

$$
A = I+\frac{i\Delta t}{2}H,
\qquad
B = I-\frac{i\Delta t}{2}H,
$$

then each update solves

$$
A\psi^{n+1}=B\psi^n.
$$

## 10) Dirichlet boundary treatment

For Dirichlet boundaries, endpoint values are fixed (typically zero):

$$
\psi(x_{\min},t)=\psi(x_{\max},t)=0.
$$

Implementation-wise, the solver operates on interior points and then writes them back into a full array whose boundaries remain fixed. This keeps the matrix tridiagonal for the evolved unknowns.

## 11) Periodic boundary treatment

For periodic boundaries,

$$
\psi(x_{\min},t)=\psi(x_{\max},t),
$$

so the discrete operator couples the first and last nodes.

## 12) Why periodic boundaries introduce corner couplings

In matrix form, periodic wrapping adds nonzero corner entries linking indices $0$ and $N-1$. The matrix is no longer strictly tridiagonal in the simple Dirichlet sense, because nearest-neighbor connectivity wraps across the boundary.

## 13) Dense solver workflow

The dense workflow is:

1. Assemble $A$ and $B$ for the chosen boundary condition.
2. Compute right-hand side $b=B\psi^n$.
3. Solve linear system with a generic dense solver (`numpy.linalg.solve` or equivalent path in the script).
4. Update state and diagnostics.

This is general and convenient, and it also handles periodic corner couplings directly.

## 14) Thomas algorithm workflow

For Dirichlet tridiagonal systems, the Thomas workflow is:

1. Extract subdiagonal, diagonal, and superdiagonal arrays.
2. Build RHS $b$ from the Crank–Nicolson right-hand side.
3. Apply forward elimination.
4. Apply backward substitution.
5. Insert interior solution into full state with fixed boundaries.

This reduces solve cost by exploiting tridiagonal structure.

## 15) Norm calculation

The code repeatedly computes

$$
N(t)=\int |\psi(x,t)|^2 dx
$$

as a numerical integral over the grid. The norm trace is used as a primary physical-consistency diagnostic.

## 16) Convergence checks

Convergence workflows run simulations with varying discretization (`dx`, `dt`) and compare error trends across:

- solver choice (dense vs Thomas where applicable), and
- boundary condition (Dirichlet vs periodic).

Generated figures summarize how errors change under refinement.

## 17) Computational cost checks

Performance workflows time repeated runs for different grid sizes and solver/boundary combinations, then plot computational cost curves. These figures quantify where tridiagonal structure gives practical speedups.

## 18) Analytical/spectral comparison workflow

For selected cases (e.g., infinite well, harmonic contexts), the script compares Crank–Nicolson outputs with analytical or spectral reference evolution to verify phase/amplitude behavior and absolute error trends.

## 19) Reflection/transmission estimation

In scattering scenarios, reflected and transmitted probabilities are estimated by integrating $|\psi|^2$ over left/right regions after interaction:

$$
R(t)=\int_{\text{left}} |\psi|^2 dx,
\qquad
T(t)=\int_{\text{right}} |\psi|^2 dx.
$$

These diagnostics are plotted for finite wells and barriers.

## 20) GIF and figure generation

The script uses Matplotlib plotting/animation to save:

- wave-packet evolution GIFs,
- norm/convergence/performance PNG diagnostics,
- analytical-comparison visualizations,
- scattering reflection/transmission plots.

Outputs are organized under `figures/` by topic (`free_packet`, `validation`, `performance`, `analytical_comparison`, `scattering`, etc.).

## 21) Practical implementation summary

The implementation follows a full numerical workflow:

1. build grid and potential,
2. initialize and normalize complex state,
3. assemble Crank–Nicolson operators,
4. evolve with dense or tridiagonal solve,
5. compute physical/numerical diagnostics,
6. save figures and animations for interpretation.

This keeps method, validation, and reporting tightly connected in a single reproducible script.

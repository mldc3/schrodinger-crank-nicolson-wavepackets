# Numerical Method

The simulations solve the 1D time-dependent Schrödinger equation using a finite-difference spatial grid and fixed time steps. The wave function is stored as complex values $\psi_j^n \approx \psi(x_j,t^n)$ on a uniform grid.

The spatial second derivative is approximated with the centred stencil $\partial_x^2\psi \approx (\psi_{j+1}-2\psi_j+\psi_{j-1})/\Delta x^2$. This gives a discrete Hamiltonian with nearest-neighbour coupling.

Time evolution is performed with the Crank–Nicolson method. Writing the equation as $i\partial_t\psi=H\psi$, the update is $(I+\frac{i\Delta t}{2}H)\psi^{n+1}=(I-\frac{i\Delta t}{2}H)\psi^n$, or equivalently $A\psi^{n+1}=B\psi^n$.

This method is used because it is stable and preserves the probability norm much better than a simple explicit scheme. Since quantum evolution should conserve $\int |\psi|^2dx$, norm conservation is one of the main validation checks.

The code supports two boundary conditions. Dirichlet boundaries fix the wave function to zero at the endpoints, representing hard walls. Periodic boundaries connect the first and last grid points, representing a ring-like domain.

For Dirichlet boundaries, the Crank–Nicolson matrix is tridiagonal. The code can solve this system either with a dense linear solver or with the Thomas algorithm, which is more efficient for tridiagonal systems.

For periodic boundaries, the first and last grid points are coupled, so the matrix is not purely tridiagonal in the simple form used here. The current implementation therefore uses dense linear algebra for the periodic case.

The repository includes diagnostics for norm conservation, convergence with respect to grid and timestep size, computational-cost comparisons between solvers, analytical or spectral validation cases, and reflection/transmission analysis for scattering potentials.

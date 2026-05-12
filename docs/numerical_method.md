# Theory Background and Numerical Method: 1D Time-Dependent Schrödinger Equation

This document explains the physical background, numerical method, implementation choices, and diagnostic tests used in the one-dimensional Schrödinger wave-packet simulations.

The project studies the time evolution of a quantum wave function $\psi(x,t)$ under different potentials and boundary conditions. The main numerical method is the Crank-Nicolson scheme, which is used because it is stable and preserves the quantum norm very well. The code also compares different linear solvers, mainly a dense solver and the Thomas algorithm for tridiagonal systems.

The general computational workflow is:

```python
# 1. Build spatial grid
x = np.linspace(x_min, x_max, N)
dx = x[1] - x[0]

# 2. Choose time step and number of steps
dt = ...
steps = ...

# 3. Define potential V(x)
V = V_func(x)

# 4. Build initial Gaussian wave packet
psi0 = gaussian_packet(x, x0, sigma, k0)

# 5. Normalize wave function
psi0 = psi0 / sqrt(integral(|psi0|^2 dx))

# 6. Evolve with Crank-Nicolson
for n in range(steps):
    psi = CN_step(psi, V, dx, dt, resolver="thomas", cc_tipo="dirichlet")

# 7. Compute diagnostics
norm = integral(|psi|^2 dx)
density = |psi|^2
reflection, transmission = integrate_left_right_regions(...)

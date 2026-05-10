# Results Summary

This document summarizes the main numerical results of the one-dimensional Schrödinger simulations. The focus is on qualitative validation: correct wave-packet evolution, conservation of probability, comparison between solvers, convergence behaviour and physically meaningful reflection/transmission diagnostics.

---

## 1. Free Gaussian wave packet

The free-particle case uses

$$
V(x)=0.
$$

The initial state is a normalized Gaussian packet with a complex phase factor. The phase gives the packet a nonzero mean momentum, so the probability density moves across the numerical domain.

The free-packet animation shows the expected behaviour: the packet propagates while also spreading. This spreading is a purely quantum effect. A localized wave packet is not a single momentum eigenstate; it contains a distribution of momenta, and the different components accumulate different phases during time evolution.

Recommended figure:

![Free Gaussian packet](../figures/free_packet/free_gaussian_packet_dirichlet.gif)

---

## 2. Boundary-condition comparison

The repository compares Dirichlet and periodic boundary conditions.

With Dirichlet boundary conditions, the wave function vanishes at the edges. This behaves like a hard-wall box. If the wave packet reaches the boundary, it reflects.

With periodic boundary conditions, the packet leaving one side of the domain reappears from the other side. This tests a different physical topology and also changes the matrix structure of the numerical problem.

Recommended figure:

![Periodic free packet](../figures/free_packet/free_gaussian_packet_periodic.gif)

The comparison is useful because it shows that boundary conditions are not just a numerical detail. They define the physical problem being solved.

---

## 3. Dense solver versus Thomas algorithm

For Dirichlet boundary conditions, the Crank–Nicolson matrices are tridiagonal. This allows the linear system to be solved with the Thomas algorithm.

The dense solver and the Thomas solver should produce essentially the same physical evolution when applied to the same Dirichlet problem. The purpose of comparing them is not to change the physics, but to test whether the specialized tridiagonal solver reproduces the same result more efficiently.

Recommended figure:

![CN dense versus Thomas](../figures/boundary_conditions/free_packet_dirichlet_cn_dense_vs_thomas.gif)

The important conclusion is that exploiting matrix structure matters. In a one-dimensional finite-difference problem, the Hamiltonian couples only nearest neighbours, so a tridiagonal solver is mathematically natural and computationally efficient.

---

## 4. Conservation of norm

The norm is

$$
N(t)=\int |\psi(x,t)|^2 dx.
$$

For Schrödinger evolution, this quantity should remain constant. The norm-conservation plot is therefore one of the most important validation tests.

Recommended figure:

![Norm conservation](../figures/validation/norm_conservation_free_packet.png)

The expected behaviour is that $N(t)$ remains very close to its initial value. Small deviations can arise from finite grid resolution, floating-point roundoff, boundary effects or explicit renormalization choices in some parts of the code.

A method that produces a large systematic norm drift would not be reliable for quantum dynamics. Crank–Nicolson is chosen precisely because it has strong norm-conservation properties for Hermitian Hamiltonians.

---

## 5. Convergence diagnostics

The convergence plots compare how the numerical error changes with the spatial step, time step, boundary condition and solver choice.

Recommended figures:

![Dense periodic convergence](../figures/validation/convergence_cn_dense_periodic.png)

![Thomas Dirichlet convergence](../figures/validation/convergence_cn_thomas_dirichlet.png)

![Dense Dirichlet convergence](../figures/validation/convergence_cn_dense_dirichlet.png)

The main interpretation is that numerical accuracy depends on both $\Delta x$ and $\Delta t$. Refining the grid generally improves the representation of the second derivative, while reducing the time step improves the temporal approximation.

Crank–Nicolson is second order in time, but the full observed error also depends on spatial discretization, boundary conditions and the diagnostic used to measure error.

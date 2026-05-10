# Numerical Method

This document describes the numerical implementation used in the repository.

---

## 1. Spatial grid

The simulation uses a one-dimensional uniform grid:

$$
x_j=x_{\min}+j\Delta x.
$$

The wave function is stored as a complex NumPy array:

$$
\psi_j^n \approx \psi(x_j,t^n).
$$

The potential is also stored on the same grid:

$$
V_j=V(x_j).
$$

---

## 2. Initial condition

The default initial state is a Gaussian wave packet:

$$
\psi(x,0)
=
\exp\left[-\frac{(x-x_0)^2}{2\sigma^2}\right]
\exp(ik_0x).
$$

The code normalizes the state using numerical integration so that

$$
\int |\psi(x,0)|^2dx=1.
$$

---

## 3. Crank–Nicolson update

The discrete Schrödinger equation is evolved using Crank–Nicolson:

$$
A\psi^{n+1}=B\psi^n.
$$

The matrices $A$ and $B$ are built from the finite-difference Hamiltonian. The method is implicit, so each time step requires solving a linear system.

---

## 4. Dirichlet boundary conditions

For Dirichlet boundaries, the endpoints satisfy

$$
\psi_0=\psi_{N-1}=0.
$$

Only the interior points are evolved. This produces a tridiagonal linear system. The code can solve it either with a dense solver or with the Thomas algorithm.

---

## 5. Periodic boundary conditions

For periodic boundaries, the first and last grid points are coupled. This creates additional corner entries in the matrix. The system is no longer strictly tridiagonal, so the simple Thomas implementation is not used for this case.

---

## 6. Thomas algorithm

The Thomas algorithm solves tridiagonal systems in linear time. It is used when the Crank–Nicolson matrix is tridiagonal, especially for Dirichlet boundary conditions.

The algorithm has two stages:

1. forward elimination,
2. backward substitution.

This makes it much more efficient than dense solvers for large one-dimensional grids.

---

## 7. Norm diagnostic

At each time, the norm is

$$
N(t)=\int |\psi(x,t)|^2dx.
$$

This is computed numerically. Since Schrödinger evolution conserves probability, the norm should remain close to one.

---

## 8. Reference solutions

Some simulations are compared with analytical or spectral reference solutions.

For the infinite well, analytical eigenstates are used. For more general potentials, the Hamiltonian can be diagonalized numerically to build a spectral reference solution:

$$
\psi(x,t)=\sum_n c_n\phi_n(x)e^{-iE_nt}.
$$

---

## 9. Scattering diagnostics

For barrier and finite-well scattering, the code estimates reflected and transmitted probability by integrating $|\psi|^2$ on opposite sides of the potential region.

This gives qualitative reflection and transmission diagnostics for the evolved packet.

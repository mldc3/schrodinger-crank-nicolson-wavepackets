# Theory

This project studies the one-dimensional time-dependent Schrödinger equation:

\[
i\,\frac{\partial \psi(x,t)}{\partial t} = \left(-\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x)\right)\psi(x,t).
\]

The wave function \(\psi\) is complex-valued, and the probability density is \(|\psi|^2\). The total probability (norm) is

\[
\int |\psi(x,t)|^2\,dx,
\]

which should remain constant under unitary evolution.

The repository explores physically relevant settings present in the source implementation:
- Free-packet propagation \((V=0)\),
- Infinite-well and hard-wall behavior,
- Harmonic-like confinement,
- Finite wells and barriers for scattering diagnostics.

Both Dirichlet and periodic boundary conditions are included in the numerical experiments.

---

## 7. Crank–Nicolson scheme

The Crank–Nicolson method is obtained by averaging the Hamiltonian action between the old and new time levels. If the Schrödinger equation is written as

$$
i\frac{\partial \psi}{\partial t}=H\psi,
$$

then Crank–Nicolson gives

$$
\left(I+\frac{i\Delta t}{2}H\right)\psi^{n+1}
=
\left(I-\frac{i\Delta t}{2}H\right)\psi^n.
$$

This can be written compactly as

$$
A\psi^{n+1}=B\psi^n.
$$

At each time step, the method solves a linear system. This is more expensive than an explicit update, but it has much better stability and norm-conservation properties.

For a Hermitian Hamiltonian, Crank–Nicolson is unitary up to numerical roundoff. This makes it especially appropriate for quantum mechanics.

---

## 8. Matrix structure

Using centred finite differences, the discrete Hamiltonian contains nearest-neighbour couplings:

$$
H_{j,j}=\frac{2}{\Delta x^2}+V_j,
$$

$$
H_{j,j+1}=H_{j,j-1}=-\frac{1}{\Delta x^2}.
$$

For Dirichlet boundary conditions, the interior problem is tridiagonal. This is important because a tridiagonal system can be solved efficiently using the Thomas algorithm.

The matrices $A$ and $B$ inherit this structure. The code therefore compares a generic dense solver with a specialized tridiagonal solver.

---

## 9. Thomas algorithm

The Thomas algorithm is a specialized direct solver for tridiagonal systems. A tridiagonal matrix contains nonzero entries only on the main diagonal, the upper diagonal and the lower diagonal.

The algorithm performs:

1. forward elimination,
2. backward substitution.

Its computational cost scales as

$$
O(N),
$$

whereas dense linear solvers generally scale much worse for large matrices.

This is why Thomas is attractive for one-dimensional finite-difference quantum simulations with Dirichlet boundary conditions.

---

## 10. Boundary conditions

Two types of boundary conditions are considered.

Dirichlet boundary conditions impose

$$
\psi(x_{\min},t)=0,
\qquad
\psi(x_{\max},t)=0.
$$

These represent hard walls. They are appropriate for an infinite well or for a finite numerical box where the wave packet is not expected to reach the boundaries too strongly.

Periodic boundary conditions impose

$$
\psi(x_{\min},t)=\psi(x_{\max},t).
$$

They represent a ring-like topology. The matrix is no longer purely tridiagonal because the first and last grid points are coupled. For that reason, the Thomas algorithm is not directly applicable to the periodic case in the same simple form.

---

## 11. Norm conservation

The norm is computed as

$$
N(t)=\int |\psi(x,t)|^2dx.
$$

Numerically, this is approximated by a quadrature or grid sum. Since $|\psi|^2$ is a probability density, the norm should remain close to one.

Norm conservation is one of the most important validation checks in this project. If the norm drifts significantly, possible causes include:

- unstable time integration,
- excessively large time step,
- inconsistent boundary treatment,
- loss of unitarity from the numerical method,
- accumulated roundoff error.

Crank–Nicolson should preserve the norm very accurately when implemented correctly.

---

## 12. Analytical and spectral reference solutions

Numerical results are compared with reference solutions when possible.

For an infinite well, stationary eigenstates have sinusoidal spatial dependence and phase evolution

$$
\psi_n(x,t)=\phi_n(x)e^{-iE_nt}.
$$

A general initial state can be expanded as

$$
\psi(x,0)=\sum_n c_n\phi_n(x).
$$

The time evolution is then

$$
\psi(x,t)=\sum_n c_n\phi_n(x)e^{-iE_nt}.
$$

For generic finite wells or barriers, a discrete Hamiltonian can be diagonalized to obtain a spectral reference solution. This gives a useful independent check of the Crank–Nicolson evolution.

---

## 13. Reflection and transmission

For scattering problems, the packet approaches a finite well or barrier. Part of the probability density may be reflected and part may be transmitted.

The reflected probability can be estimated by integrating on the left side:

$$
R(t)=\int_{x

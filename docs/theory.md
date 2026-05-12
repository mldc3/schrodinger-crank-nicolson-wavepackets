# Theory Background: 1D Time-Dependent Schrödinger Equation

This document explains the physical and numerical background of the one-dimensional Schrödinger simulations implemented in this repository.

The project studies the time evolution of quantum wave packets under different potentials and boundary conditions. The main numerical method is the Crank–Nicolson scheme, combined with either a dense linear solver or the Thomas algorithm when the matrix structure is tridiagonal.

---

## 1. Physical problem

The time-dependent Schrödinger equation describes the evolution of a quantum state. In one spatial dimension, and in the dimensionless units used in this project, it can be written as

$$
i\frac{\partial \psi(x,t)}{\partial t} = -\frac{\partial^2 \psi(x,t)}{\partial x^2}
+
V(x)\psi(x,t).
$$

Here, $\psi(x,t)$ is the complex wave function and $V(x)$ is the external potential. The measurable probability density is

$$
\rho(x,t)=|\psi(x,t)|^2.
$$

The total probability must remain normalized:

$$
\int |\psi(x,t)|^2 dx = 1.
$$

A good numerical method for Schrödinger evolution should therefore preserve the norm as accurately as possible.

---

## 2. Why this equation is different from diffusion

Although the Schrödinger equation contains a second spatial derivative, it is not a diffusion equation. The factor $i$ in front of the time derivative changes the physical character of the equation.

A diffusion equation smooths profiles and dissipates information. The Schrödinger equation instead generates unitary time evolution. Probability density can spread, interfere, reflect and transmit, but the total probability should remain conserved.

This is why methods that work for the heat equation are not automatically appropriate for quantum evolution. Numerical stability and norm conservation become central diagnostics.

---

## 3. Wave packets

A common initial state is a Gaussian wave packet:

$$
\psi(x,0) =
\exp\left[-\frac{(x-x_0)^2}{2\sigma^2}\right]
\exp(ik_0x).
$$

The Gaussian envelope localizes the particle around $x_0$, while the complex phase gives it a mean momentum related to $k_0$.

After constructing the packet, it is normalized numerically:

$$
\psi(x,0)
\leftarrow
\frac{\psi(x,0)}
{\sqrt{\int |\psi(x,0)|^2 dx}}.
$$

This makes the total probability equal to one at the start of the simulation.

---

## 4. Potentials studied

The numerical code can evolve the wave packet under several representative potentials.

For a free particle,

$$
V(x)=0.
$$

For an infinite well, the wave function is confined by hard-wall boundary conditions. In the interior, the potential may be treated as zero, while the walls are imposed through Dirichlet conditions.

For a harmonic oscillator,

$$
V(x)=\frac{1}{2}\omega^2x^2.
$$

For a finite well,

$$
V(x)=
\begin{cases}
-V_0, & |x|<a/2,\\
0, & |x|\ge a/2.
\end{cases}
$$

For a rectangular barrier,

$$
V(x)=
\begin{cases}
V_0, & |x|<a/2,\\
0, & |x|\ge a/2.
\end{cases}
$$

These cases test qualitatively different quantum behaviour: free spreading, confinement, oscillation, reflection and transmission.

---

## 5. Spatial and temporal discretization

The spatial domain is discretized as

$$
x_j=x_{\min}+j\Delta x,
$$

and time is discretized as

$$
t^n=n\Delta t.
$$

The wave function on the grid is written as

$$
\psi_j^n \approx \psi(x_j,t^n).
$$

The second spatial derivative is approximated by the centred finite difference

$$
\frac{\partial^2 \psi}{\partial x^2}
\approx
\frac{
\psi_{j+1}-2\psi_j+\psi_{j-1}
}{\Delta x^2}.
$$

This approximation naturally produces a tridiagonal matrix for Dirichlet boundary conditions.

---

## 6. Explicit FTCS and its limitation

A direct explicit discretization updates the wave function using only the previous time level. This is simple, but it is not the preferred method for the Schrödinger equation.

The reason is that the exact quantum evolution is unitary. A poor explicit scheme can introduce artificial growth or damping of the wave-function norm. Even if the equation is physically conservative, the numerical method may not be.

Therefore, the project focuses on Crank–Nicolson, which is much better suited to Schrödinger evolution.

---

## 7. Crank–Nicolson scheme

The Crank–Nicolson method is obtained by averaging the Hamiltonian action between the old and new time levels. If the Schrödinger equation is written as

$$
i\frac{\partial \psi}{\partial t}=H\psi,
$$

then Crank–Nicolson gives

$$
\left(I+\frac{i\Delta t}{2}H\right)\psi^{n+1} =
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
R(t)=\int_{x<x_L}|\psi(x,t)|^2dx.
$$

The transmitted probability can be estimated by integrating on the right side:

$$
T(t)=\int_{x>x_R}|\psi(x,t)|^2dx.
$$

In an ideal large domain after the scattering event, one expects approximately

$$
R+T \approx 1.
$$

This provides a physical diagnostic of probability conservation and scattering behaviour.

---

## 14. Summary

The project combines several central ideas from computational quantum mechanics:

- finite-difference discretization of the Schrödinger equation,
- Gaussian wave-packet initialization,
- Crank–Nicolson time evolution,
- tridiagonal linear systems,
- Thomas algorithm acceleration,
- Dirichlet and periodic boundary conditions,
- norm conservation,
- comparison with analytical or spectral solutions,
- wave-packet spreading,
- reflection and transmission.

Together, these elements form a compact but complete numerical quantum-dynamics project.

# Theory Background: 1D Time-Dependent Schrödinger Equation

This document explains the physical and numerical background of the one-dimensional Schrödinger simulations implemented in this repository.

The project studies the time evolution of quantum wave packets under different potentials and boundary conditions. The main numerical method is the Crank–Nicolson scheme, combined with either a dense linear solver or the Thomas algorithm when the matrix structure is tridiagonal.

---

## 1. Physical problem

The time-dependent Schrödinger equation describes the evolution of a quantum state. In one spatial dimension, and in the dimensionless units used in this project, it can be written as

$$
i\frac{\partial \psi(x,t)}{\partial t}
=
-\frac{\partial^2 \psi(x,t)}{\partial x^2}
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
\psi(x,0)
=
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

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

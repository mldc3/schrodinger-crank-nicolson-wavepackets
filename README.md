# 1D Time-Dependent Schrödinger Equation with Crank–Nicolson Wave Packets

This repository is a cleaned, portfolio-ready version of a computational physics project for the 1D time-dependent Schrödinger equation. It studies Gaussian wave-packet dynamics with finite differences and Crank–Nicolson time stepping across multiple physical scenarios.

## What equation is solved

The simulations evolve the complex wave function $\psi(x,t)$ using

$$
i\frac{\partial \psi}{\partial t} = -\frac{\partial^2 \psi}{\partial x^2} + V(x)\psi.
$$

## Why this equation matters

The time-dependent Schrödinger equation is the core dynamical model of nonrelativistic quantum mechanics. It predicts how probability amplitudes propagate, interfere, reflect, and transmit in external potentials.

## What Crank–Nicolson contributes

Crank–Nicolson provides an implicit, stable, second-order time update with strong norm-conservation behavior for Hermitian Hamiltonians. In practice, this makes long-time quantum propagation more reliable than naive explicit updates.

## Boundary conditions studied

- **Dirichlet boundaries:** hard-wall endpoints, useful for boxed domains and infinite-well-style setups.
- **Periodic boundaries:** wrap-around domain, where left and right edges are coupled.

## Why norm conservation matters

A central diagnostic is the discrete probability norm

$$
N(t)=\int |\psi(x,t)|^2\,dx.
$$

For physically consistent Schrödinger evolution, $N(t)$ should remain approximately constant. The repository includes dedicated norm-conservation plots.

![Norm conservation (free packet)](figures/validation/norm_conservation_free_packet.png)

## What the Thomas algorithm comparison demonstrates

For Dirichlet finite-difference systems, the matrix is tridiagonal. The project compares:

- dense linear solves, and
- Thomas tridiagonal solves.

This demonstrates that exploiting matrix structure can preserve the same physics with lower computational cost.

## What the analytical comparisons validate

The infinite-well and harmonic-oscillator comparisons test whether numerical propagation reproduces known analytical/spectral behavior and expected error trends.

![Infinite well: CN vs analytical](figures/analytical_comparison/infinite_well_cn_vs_analytical.gif)
![Harmonic packet: CN vs analytical](figures/analytical_comparison/harmonic_packet_cn_vs_analytical.gif)

## What the scattering figures show

Finite-well and barrier cases show wave-packet splitting into reflected and transmitted components, with diagnostics based on integrated probability on each side of the interaction region.

![Rectangular barrier reflection/transmission](figures/scattering/rectangular_barrier_reflection_transmission.png)
![Finite well reflection/transmission](figures/scattering/finite_well_reflection_transmission.png)

## How to run

1. Create and activate a Python environment.
2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Run the main script:

```bash
python src/schrodinger_crank_nicolson.py
```

## Where things are

- Main code: `src/schrodinger_crank_nicolson.py`
- Documentation: `docs/theory.md`, `docs/numerical_method.md`, `docs/results_summary.md`, `docs/sources_and_notes.md`
- Figures: `figures/free_packet/`, `figures/boundary_conditions/`, `figures/validation/`, `figures/analytical_comparison/`, `figures/scattering/`, `figures/performance/`

## Skills demonstrated

- Finite-difference PDE discretization for quantum dynamics
- Implicit time integration (Crank–Nicolson)
- Complex-valued linear algebra and tridiagonal solvers
- Boundary-condition analysis (Dirichlet vs periodic)
- Numerical validation (norm, convergence, analytical checks)
- Scientific communication with organized results and reproducible documentation

## Author

**María Lourdes Domínguez Cacho**  
Final-semester Physics student, University of Alicante  
GitHub: [mldc3](https://github.com/mldc3)

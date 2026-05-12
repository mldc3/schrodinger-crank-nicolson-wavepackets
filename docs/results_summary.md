# Results Summary and Physical Interpretation

## 1. General interpretation of the numerical results

The simulations show that the Crank–Nicolson method is well suited for one-dimensional quantum time evolution. The method produces stable wave-packet propagation, preserves the norm well, and allows comparisons between different boundary conditions and solvers.

The main physical behaviour observed in the figures is consistent with quantum mechanics:

- a free Gaussian wave packet moves because of its initial momentum,
- the same packet spreads because it contains a distribution of momentum components,
- Dirichlet boundaries produce reflection from hard walls,
- periodic boundaries allow the packet to re-enter from the other side,
- bound systems such as the infinite well and harmonic oscillator show confined motion,
- barriers and wells split the packet into reflected and transmitted parts.

The main numerical behaviour is also consistent with the theory:

- Crank–Nicolson gives stable time evolution,
- the Thomas algorithm reproduces the dense solver for Dirichlet tridiagonal systems,
- convergence improves when the spatial and temporal grids are refined,
- dense solvers are more general but computationally more expensive,
- exploiting tridiagonal structure is important for efficiency.

---

## 2. Free Gaussian wave packet

The free-particle simulation uses the potential $V(x)=0$.

The initial state is a normalized Gaussian wave packet with a complex phase factor. It has the general form $\psi(x,0)=\exp[-(x-x_0)^2/(2\sigma^2)]\exp(ik_0x)$. The Gaussian part localizes the particle around $x_0$, while the complex phase $\exp(ik_0x)$ gives the packet a nonzero average momentum.

Recommended figure:

![Free Gaussian packet](../figures/free_packet/free_gaussian_packet_dirichlet.gif)

The animation shows the expected free-particle behaviour. The packet moves across the numerical domain, showing that the phase factor has correctly introduced momentum. At the same time, the packet spreads. This spreading is not a numerical error; it is a physical consequence of the Schrödinger equation.

A localized Gaussian wave packet is not a single momentum eigenstate. Instead, it is a superposition of different momentum components. Each component evolves with a different phase in time, so the shape of the packet changes. This produces dispersion, meaning that the width of the packet increases as time passes.

This result is an important first validation because it shows that the code is not only transporting the packet but also reproducing the dispersive nature of quantum evolution.

---

## 3. Packet spreading

The spreading of the free Gaussian packet is analysed more directly by tracking its width.

Recommended figure:

![Gaussian packet spreading](../figures/free_packet/gaussian_packet_spreading_width.png)

The width of the packet increases with time, as expected for a free quantum particle. This confirms that the simulation captures a key property of the time-dependent Schrödinger equation: even in the absence of forces, a localized wave packet does not remain rigid.

Physically, this happens because localization in position implies uncertainty in momentum. The different momentum components travel with different phase velocities and group behaviour, causing the wave packet to broaden.

Numerically, this is also a useful check. If the packet only translated without changing shape, the simulation would not be resolving the quantum dispersion correctly. If the packet spread too quickly or irregularly, it could indicate poor spatial resolution, too large a timestep, or boundary effects.

The spreading plot therefore validates both the physical model and the numerical discretization.

---

## 4. Boundary-condition comparison

The repository compares Dirichlet and periodic boundary conditions. This is important because boundary conditions are not only a technical detail; they define the physical system being simulated.

With Dirichlet boundary conditions, the wave function satisfies $\psi(x_{\min},t)=0$ and $\psi(x_{\max},t)=0$. This represents hard walls at the edges of the domain. If the wave packet reaches a boundary, it reflects.

With periodic boundary conditions, the wave function satisfies $\psi(x_{\min},t)=\psi(x_{\max},t)$. This means that the numerical domain behaves like a ring. A packet leaving one side of the domain reappears from the other side.

Recommended figure:

![Periodic free packet](../figures/free_packet/free_gaussian_packet_periodic.gif)

The difference between the two cases is physically meaningful. In the Dirichlet case, the simulation resembles a particle in a box with walls. In the periodic case, the simulation resembles motion on a closed loop. Therefore, even with the same initial wave packet and potential, the long-time behaviour is different.

This comparison also affects the matrix structure. Dirichlet boundary conditions lead to a tridiagonal interior system, while periodic boundary conditions introduce corner couplings between the first and last grid points. This is why the Thomas algorithm is directly suitable for Dirichlet systems but not for the periodic case in its basic form.

The boundary-condition comparison confirms that the code correctly implements two different physical settings.

---

## 5. Dense solver versus Thomas algorithm

For Dirichlet boundary conditions, the finite-difference Hamiltonian only couples neighbouring grid points. This produces a tridiagonal matrix in the Crank–Nicolson linear system.

The dense solver and Thomas solver are compared for the same Dirichlet problem.

Recommended figure:

![CN dense versus Thomas](../figures/boundary_conditions/free_packet_dirichlet_cn_dense_vs_thomas.gif)

The two curves should overlap almost perfectly. This is expected because both solvers solve the same mathematical system. The physical evolution should not depend on whether the linear system is solved by a general dense solver or by a tridiagonal algorithm.

The purpose of this comparison is therefore numerical validation. If the Thomas result differs significantly from the dense result, it would indicate an error in the extraction of the diagonals, the construction of the matrices, or the Thomas implementation.

The result shows that the Thomas algorithm reproduces the dense solver while being more appropriate for tridiagonal systems. This is important because the Thomas algorithm has computational cost $O(N)$, whereas dense solving is much more expensive for large matrices.

This comparison demonstrates that using the structure of the problem matters. In one-dimensional finite differences with Dirichlet boundaries, the Hamiltonian has nearest-neighbour coupling, so a tridiagonal solver is the natural choice.

---

## 6. Conservation of norm

The probability norm is defined as $N(t)=\int |\psi(x,t)|^2 dx$.

For exact Schrödinger evolution, the norm must remain constant. This is because the time evolution is unitary: probability is not created or destroyed.

Recommended figure:

![Norm conservation](../figures/validation/norm_conservation_free_packet.png)

The norm-conservation plot is one of the most important validation tests in the project. The expected behaviour is that $N(t)$ remains very close to its initial value throughout the simulation.

Small deviations can occur because of:

- finite grid spacing,
- floating-point roundoff,
- numerical quadrature error,
- boundary treatment,
- accumulated linear-solver error,
- explicit renormalization in some comparison sections.

However, these deviations should remain small. A large monotonic drift in the norm would indicate that the method is not correctly preserving probability.

Crank–Nicolson is chosen precisely because it has strong norm-conservation properties. For a Hermitian Hamiltonian, the Crank–Nicolson update is unitary up to numerical precision. This makes it much more suitable for the Schrödinger equation than a naive explicit method.

The norm plot therefore confirms that the implementation is reliable for quantum time evolution.

---

## 7. Convergence diagnostics

The convergence plots study how the numerical error depends on the spatial step $\Delta x$, the time step $\Delta t$, the solver, and the boundary condition.

Recommended figures:

![Dense periodic convergence](../figures/validation/convergence_cn_dense_periodic.png)

![Thomas Dirichlet convergence](../figures/validation/convergence_cn_thomas_dirichlet.png)

![Dense Dirichlet convergence](../figures/validation/convergence_cn_dense_dirichlet.png)

The main idea is that a numerical method should become more accurate when the discretization is refined. Reducing $\Delta x$ improves the approximation of the second derivative, while reducing $\Delta t$ improves the time integration.

Crank–Nicolson is second order in time, and the centred finite-difference Laplacian is second order in space. Therefore, the general expectation is that the error should decrease as the grid is refined.

In practice, the observed convergence depends on the diagnostic used. In this repository, one important diagnostic is the error in norm conservation. Because Crank–Nicolson already conserves the norm very well, the errors may be very small and can be influenced by roundoff, solver precision and boundary effects.

The convergence plots are still useful because they show whether the method behaves consistently. If the error increased when using a smaller timestep or finer grid, that would indicate an implementation issue. The expected behaviour is that refined discretizations produce equal or better accuracy.

---

## 8. Computational cost comparison

The computational-cost plots compare the runtime of the dense solver and the Thomas algorithm for different discretizations.

Recommended figures:

![Dense periodic cost](../figures/performance/computational_cost_cn_dense_periodic.png)

![Thomas Dirichlet cost](../figures/performance/computational_cost_cn_thomas_dirichlet.png)

![Dense Dirichlet cost](../figures/performance/computational_cost_cn_dense_dirichlet.png)

The dense solver is general and easy to use, but it does not exploit the special structure of the matrix. It treats the matrix as a full system, even though most entries are zero.

The Thomas algorithm is specialized for tridiagonal matrices. For a one-dimensional finite-difference Hamiltonian with Dirichlet boundaries, only the main diagonal, upper diagonal and lower diagonal are nonzero. Therefore, Thomas solves the problem much more efficiently.

The expected computational behaviour is:

- dense solving becomes expensive as the number of grid points increases,
- Thomas scales much better for Dirichlet tridiagonal systems,
- periodic boundary conditions require dense solving in the current implementation because endpoint coupling breaks the simple tridiagonal form.

The performance comparison shows why mathematical structure matters in scientific computing. Two methods can give the same physical answer, but their computational efficiency can be very different.

This is one of the main numerical lessons of the project: the best method is not only the one that is accurate, but also the one that uses the structure of the discretized equations efficiently.

---

## 9. Infinite well comparison

The infinite square well is a strong validation problem because analytical solutions are known.

In an infinite well, the wave function must vanish at the boundaries. The stationary states are sinusoidal eigenfunctions, and each eigenstate evolves only by a phase factor in time.

Recommended figures:

![Infinite well CN versus analytical](../figures/analytical_comparison/infinite_well_cn_vs_analytical.gif)

![Infinite well packet with walls](../figures/analytical_comparison/infinite_well_packet_with_walls.gif)

![Infinite well absolute error](../figures/analytical_comparison/infinite_well_absolute_error.png)

If the initial state is a single eigenstate, the probability density should remain stationary. The complex wave function changes phase, but $|\psi(x,t)|^2$ does not change shape.

If the initial state is a wave packet or a superposition of eigenstates, the probability density evolves because each eigenstate accumulates a different phase. This can produce oscillations, interference patterns and reflections from the walls.

The comparison between Crank–Nicolson and the analytical solution tests whether the numerical method reproduces the known quantum dynamics. The absolute error plot gives a direct measure of the difference between the numerical and reference solutions.

The expected result is that the numerical and analytical curves remain close, with small errors caused by finite spatial resolution, finite timestep and boundary discretization.

This is a stronger test than just checking that the animation “looks right,” because it compares the simulation against a known solution.

---

## 10. Harmonic oscillator packet

The harmonic oscillator potential is $V(x)=\frac{1}{2}\omega^2x^2$.

Recommended figure:

![Harmonic packet](../figures/analytical_comparison/harmonic_packet_cn_vs_analytical.gif)

This test is useful because the harmonic oscillator has a very clear physical interpretation. The potential confines the packet and pulls it back toward the centre. A displaced Gaussian packet should oscillate in the quadratic potential.

The harmonic oscillator is different from the free-particle case because the packet is not simply moving across the domain and spreading freely. Instead, the potential continuously changes the momentum of the packet and causes bounded motion.

A good numerical simulation should show the packet moving back and forth in the potential well. The probability density should remain localized if the state is close to a coherent-state-like Gaussian packet.

This case validates that the code can handle nonzero smooth potentials, not only the free equation. It also tests whether the potential term has been included correctly in the Crank–Nicolson matrices.

---

## 11. Analytical and spectral comparisons

The project uses analytical or spectral reference solutions to validate the numerical method.

For the infinite well, the reference solution can be constructed from known eigenfunctions.

For more general potentials such as finite wells or barriers, the code can build a discrete Hamiltonian and diagonalize it. The wave function is then expanded in the eigenvectors of this Hamiltonian.

The general idea is $\psi(x,t)=\sum_n c_n\phi_n(x)e^{-iE_nt}$, where the coefficients are obtained from the initial state as $c_n=\langle \phi_n|\psi(0)\rangle$.

This provides an independent comparison against the time-stepping Crank–Nicolson method.

These comparisons are important because they check more than norm conservation. A method could conserve the norm but still evolve with the wrong phase or wrong dispersion relation. Comparing with reference solutions checks whether the shape, phase behaviour and probability density evolve correctly.

One point to handle carefully is the Hamiltonian convention. If one part of the code uses $H=-\partial_x^2+V$ and another uses $H=-\frac{1}{2}\partial_x^2+V$, the two results will not match exactly because the kinetic term is different. Therefore, reference comparisons must use the same Hamiltonian convention as the Crank–Nicolson evolution.

---

## 12. Reflection and transmission

The finite well and rectangular barrier simulations study scattering.

In a scattering problem, a wave packet approaches a localized potential. Part of the probability density can be reflected, and part can be transmitted.

Recommended figures:

![Rectangular barrier reflection and transmission](../figures/scattering/rectangular_barrier_reflection_transmission.png)

![Finite well reflection and transmission](../figures/scattering/finite_well_reflection_transmission.png)

The reflected and transmitted probabilities are estimated by integrating the probability density on different sides of the scattering region.

For example, $R(t)=\int_{\mathrm{left}}|\psi(x,t)|^2dx$ and $T(t)=\int_{\mathrm{right}}|\psi(x,t)|^2dx$.

These quantities are useful diagnostics, but they must be interpreted carefully. At early times, the wave packet may still be located on the left side before it has interacted with the barrier or well. Therefore, the left-side probability is not initially “reflected probability” in the physical scattering sense. It becomes meaningful as reflected probability only after the packet has interacted with the potential and the reflected and transmitted parts have separated.

The same applies to transmission. The transmitted probability becomes meaningful after the right-moving part has passed through the scattering region.

A good scattering simulation should show:

- initially, most probability located on the incident side,
- during interaction, probability density overlaps with the potential region,
- after interaction, the packet separates into reflected and transmitted components,
- the sum of probabilities remains approximately conserved.

In an ideal sufficiently large domain after scattering, one expects $R+T\approx 1$, up to probability still inside the interaction region and numerical errors.

---

## 13. Rectangular barrier behaviour

For a rectangular barrier, the potential is positive in a finite region. Classically, a particle with energy below the barrier would be fully reflected. Quantum mechanically, part of the wave packet can tunnel through.

Therefore, the rectangular barrier test checks whether the simulation captures a key quantum effect: partial transmission through a classically forbidden or partially forbidden region.

The amount of reflection and transmission depends on:

- barrier height,
- barrier width,
- wave-packet momentum,
- wave-packet energy spread,
- domain size,
- final simulation time.

A higher or wider barrier generally increases reflection and decreases transmission. A packet with larger mean momentum generally transmits more easily.

The plotted reflection and transmission curves should show the splitting of probability as the packet interacts with the barrier. This provides a physical validation of the potential implementation and the probability integration diagnostics.

---

## 14. Finite well behaviour

For a finite well, the potential is negative in a finite region. Instead of repelling the packet, the well attracts it.

A finite well can produce reflection, transmission and temporary trapping of probability density inside the well. This is different from the barrier case because the potential region can support bound or quasi-bound behaviour depending on the parameters.

When the packet reaches the well, part of it may transmit through, part may reflect, and part may remain temporarily localized near the well before escaping. This can make the reflection/transmission curves more complicated than for a simple barrier.

Therefore, for a finite well, $R(t)$ and $T(t)$ should be interpreted together with the probability remaining near the well region. If $R+T$ is less than 1 at intermediate times, this does not necessarily mean probability is lost. It may mean that some probability density is still inside or near the well.

This result is useful because it shows that the code can handle attractive potentials and not only repulsive barriers.

---

## 15. Interpretation of probability conservation in scattering

For scattering simulations, probability conservation should be checked globally.

The total probability is $N(t)=\int_{-\infty}^{+\infty}|\psi(x,t)|^2dx$.

The reflected and transmitted probabilities only account for selected spatial regions. Therefore, at intermediate times, $R(t)+T(t)$ may be less than 1 because part of the wave function is still in the central region.

A more complete decomposition is $N(t)=R(t)+P_{\mathrm{middle}}(t)+T(t)$, where $P_{\mathrm{middle}}(t)$ is the probability inside the scattering region.

After a sufficiently long time, when the packet has separated into left- and right-moving components, $P_{\mathrm{middle}}(t)$ should become small. Then $R(t)+T(t)$ should approach the total norm.

This interpretation is important because otherwise early-time values of $R$ and $T$ can be misleading.

---

## 16. Validation hierarchy

The results can be understood as a sequence of validation tests.

First, the free packet checks whether the code propagates and disperses a wave function correctly.

Second, norm conservation checks whether probability is preserved.

Third, the dense-versus-Thomas comparison checks whether the specialized solver reproduces the general solver.

Fourth, boundary-condition comparisons check whether Dirichlet and periodic domains are implemented correctly.

Fifth, convergence studies check whether the numerical error behaves consistently as $\Delta x$ and $\Delta t$ are changed.

Sixth, analytical and spectral comparisons check whether the numerical solution follows known or independently computed quantum dynamics.

Finally, scattering simulations check whether the code produces physically meaningful reflection and transmission.

Together, these tests provide a coherent validation workflow.

---

## 17. Main limitations of the results

The results are physically meaningful, but they should be interpreted with some limitations.

First, the simulations are one-dimensional. Real quantum systems can have higher-dimensional effects that are not represented here.

Second, the finite numerical domain can affect long-time behaviour. If a packet reaches a Dirichlet boundary, it reflects from the wall. This may interfere with the physical scattering process if the domain is not large enough.

Third, the periodic case represents a ring-like domain, not an infinite line. This is useful for testing but has a different physical interpretation.

Fourth, reflection and transmission probabilities require enough time for the reflected and transmitted packets to separate. Before that, the left/right integrals are only regional probabilities.

Fifth, analytical comparisons must use the same Hamiltonian convention as the numerical solver. A mismatch in the kinetic prefactor changes the dynamics.

Sixth, dense periodic solving is convenient but not optimal for large systems. A cyclic tridiagonal solver or sparse method would be more efficient for periodic boundary conditions.

---

## 18. Main conclusions

The results show that the Crank–Nicolson method provides a robust and stable framework for one-dimensional Schrödinger simulations.

The free-packet simulations reproduce the expected propagation and quantum spreading. This confirms that the code captures the dispersive nature of the Schrödinger equation.

The boundary-condition comparison shows that Dirichlet and periodic boundaries correspond to different physical problems. Dirichlet boundaries behave like hard walls, while periodic boundaries behave like a ring.

The dense solver and Thomas algorithm give the same physical evolution for Dirichlet problems, validating the Thomas implementation. The Thomas algorithm is preferable for large tridiagonal systems because it is much more efficient.

The norm-conservation results confirm that Crank–Nicolson preserves probability well, which is essential for quantum mechanics.

The convergence and cost studies show the expected trade-off between accuracy and computational expense. Smaller $\Delta x$ and $\Delta t$ improve the discretization, while specialized solvers reduce computational cost.

The infinite-well and harmonic-oscillator tests show that the method can reproduce known bound-system dynamics. The analytical and spectral comparisons provide stronger validation than visual inspection alone.

The barrier and finite-well simulations demonstrate scattering behaviour, including reflection, transmission and possible temporary localization near the potential region.

Overall, the repository presents a complete numerical workflow for one-dimensional quantum dynamics: initialization of wave packets, construction of the Crank–Nicolson matrices, implementation of boundary conditions, solver comparison, norm validation, convergence analysis, performance comparison and physical interpretation of scattering results.

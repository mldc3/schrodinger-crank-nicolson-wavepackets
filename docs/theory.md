# Theory Background and Numerical Method: 1D Time-Dependent Schrödinger Equation

This document explains the physical theory and numerical method used in the one-dimensional Schrödinger simulations of this repository. The project studies the time evolution of quantum wave packets under different potentials and boundary conditions using the Crank--Nicolson method. The implementation is mainly contained in `src/schrodinger_crank_nicolson.py`, where the most important parts are the construction of the spatial grid, the Gaussian initial state, the Crank--Nicolson step, the treatment of Dirichlet and periodic boundary conditions, the dense and Thomas linear solvers, the norm-conservation checks, the analytical or spectral comparisons, and the reflection/transmission diagnostics.

The aim of the code is not only to obtain animations of wave packets, but also to test if the numerical method reproduces the expected physics: conservation of probability, free spreading, reflection at boundaries, transmission through barriers, confinement in wells, and agreement with reference solutions when those are available.

---

## 1. Physical problem

The physical system is a single quantum particle moving in one spatial dimension under an external potential. The state of the particle is described by a complex wave function `psi(x,t)`. In the dimensionless convention used in the main Crank--Nicolson implementation, the time-dependent Schrödinger equation is written as $$i\frac{\partial \psi(x,t)}{\partial t}=-\frac{\partial^2\psi(x,t)}{\partial x^2}+V(x)\psi(x,t).$$

The term $-\partial_x^2\psi$ represents the kinetic-energy contribution, and the term $V(x)\psi$ represents the effect of the external potential. The potential changes the phase and shape of the wave function, and it can produce reflection, transmission, confinement or oscillatory motion depending on its form.

The probability density is $$\rho(x,t)=|\psi(x,t)|^2.$$ This gives the probability density for finding the particle near position $x$ at time $t$. Since the particle must be somewhere in the domain, the total probability should be conserved: $$\int |\psi(x,t)|^2\,dx=1.$$

A central requirement of the numerical simulation is therefore that the norm of the wave function remains close to one during the time evolution. If the norm grows or decays artificially, the numerical method is not correctly representing quantum dynamics.

---

## 2. Why the Schrödinger equation needs special care

The Schrödinger equation contains a second spatial derivative, so at first sight it may look similar to the diffusion or heat equation. However, the factor $i$ multiplying the time derivative makes the equation physically very different.

A diffusion equation dissipates information: peaks become smoother and the total profile spreads irreversibly. The Schrödinger equation instead gives unitary evolution. This means that the total probability is conserved, and the evolution is reversible in the ideal continuous problem.

Because of this, a numerical method for the Schrödinger equation should not behave like a dissipative heat-equation solver. It should avoid artificial damping and artificial growth. This is the main reason why a simple explicit method is not the best choice here, and why the project uses Crank--Nicolson.

Crank--Nicolson is implicit, time-centred, second-order accurate in time, and for a Hermitian Hamiltonian it is norm-preserving up to numerical roundoff. This makes it very appropriate for quantum-mechanical time evolution.

---

## 3. Wave packet initial condition

The simulations usually start from a Gaussian wave packet. This is a useful initial state because it is localized in space but also has a well-defined average momentum. In the code, the initial packet has the form $$\psi(x,0)=\exp\left[-\frac{(x-x_0)^2}{2\sigma^2}\right]\exp(ik_0x).$$

The Gaussian part localizes the packet around the initial position $x_0$. The parameter $\sigma$ controls the spatial width. A small $\sigma$ gives a narrow packet, while a large $\sigma$ gives a wider packet.

The complex phase factor $\exp(ik_0x)$ gives the packet a mean wave number $k_0$. In quantum mechanics, wave number is related to momentum, so this factor makes the packet move. If $k_0>0$, the packet moves mainly to the right; if $k_0<0$, it moves mainly to the left.

After creating the packet, the code normalizes it numerically using the grid. The normalization step is 

$$\psi(x,0)\leftarrow\frac{\psi(x,0)}{\sqrt{\int |\psi(x,0)|^2dx}}.$$ 

In the code this is done with a trapezoidal numerical integral such as `np.trapz(np.abs(psi0)**2, x)`.

This normalization is important because all later probabilities, such as reflection and transmission, only make physical sense if the initial total probability is one.

---

## 4. Potentials included in the project

The code is written so that different potentials can be tested by changing `V_func(x)` or by directly defining an array `V`. This is useful because different potentials test different physical behaviours and different numerical situations.

### 4.1 Free particle

For a free particle, the potential is 

$$V(x)=0.$$

In this case, the wave packet moves and spreads. There is no force from an external potential, but the packet still changes shape because the Schrödinger equation is dispersive. The packet contains several momentum components, and each component evolves with a different phase.

This case is useful as a basic validation test. The expected behaviour is propagation plus spreading, with conservation of norm.

### 4.2 Infinite well / hard-wall box

An infinite well confines the particle between hard walls. Numerically, this is usually represented using Dirichlet boundary conditions, 

$$\psi(x_{\min},t)=0,\qquad \psi(x_{\max},t)=0.$$

Inside the well, the potential can be treated as zero, while the walls are imposed by the boundary condition.

The infinite well is useful because analytical solutions are known. The stationary states have sinusoidal spatial shape and time-dependent phase. Therefore, this case provides a strong validation test for the Crank--Nicolson implementation.

### 4.3 Harmonic oscillator

The harmonic oscillator potential is 

$$V(x)=\frac{1}{2}\omega^2x^2.$$

This potential confines the packet near the centre. A displaced Gaussian packet oscillates in the potential, so this case tests whether the code can reproduce bound, oscillatory dynamics.

In the code, the harmonic-packet section compares numerical evolution with approximate analytical behaviour for the packet centre. This is useful conceptually, although the exact comparison requires consistent conventions for the Hamiltonian and packet parameters.

### 4.4 Finite square well

A finite well can be written as 

$$V(x)=\begin{cases}-V_0,& |x|<a/2,\\0,& |x|\ge a/2.\end{cases}$$

Here $V_0>0$ is the depth of the well and $a$ is its width. A wave packet incident on a finite well can be partially reflected and partially transmitted. The well can also support bound or quasi-bound behaviour depending on the energy of the packet and the depth of the well.

This potential is useful because it shows that even an attractive region can produce reflection due to wave interference and mismatch between regions.

### 4.5 Rectangular barrier

A rectangular barrier can be written as 

$$V(x)=\begin{cases}V_0,& |x|<a/2,\\0,& |x|\ge a/2.\end{cases}$$

Here $V_0>0$ is the barrier height. A wave packet incident on the barrier can split into reflected and transmitted parts. If the classical energy is below the barrier, quantum mechanics still allows nonzero transmission by tunnelling.

This is one of the most important physical examples in the project because it directly shows wave-packet scattering, reflection and transmission.

### 4.6 Double-well example

The code also contains an example of a quartic double-well potential of the form $$V(x)=V_0(x^2-a^2)^2.$$ This potential has two minima and a barrier between them. It is useful for studying confinement in two regions and tunnelling-like behaviour between wells. Even if it is not the main result, it shows that the code structure is general enough to accept different potential functions.

---

## 5. Spatial discretization

The continuous spatial domain is replaced by a uniform grid. If the domain goes from $x_{\min}$ to $x_{\max}$ and contains $N$ grid points, the grid is 

$$x_j=x_{\min}+j\Delta x,\qquad j=0,1,\ldots,N-1.$$ The grid spacing is $$\Delta x=x_{j+1}-x_j.$$

In the code this is produced by `x = np.linspace(-L, L, N)` and `dx = x[1] - x[0]`.

The wave function is stored as a complex array, where 

$$\psi_j^n\approx \psi(x_j,t^n).$$ The potential is also stored as an array, $$V_j=V(x_j).$$

The second derivative is approximated using the centred finite-difference stencil 

$$\frac{\partial^2\psi}{\partial x^2}\bigg|_{x_j}\approx\frac{\psi_{j+1}-2\psi_j+\psi_{j-1}}{\Delta x^2}.$$ 

This approximation is second-order accurate in space.

The important consequence is that each grid point only couples to its nearest neighbours. This creates a tridiagonal Hamiltonian matrix for Dirichlet boundary conditions.

---

## 6. Discrete Hamiltonian

The continuous Hamiltonian in the main convention is 

$$H=-\frac{\partial^2}{\partial x^2}+V(x).$$

After finite-difference discretization, the Hamiltonian becomes a matrix acting on the vector of wave-function values.

For an interior point, the kinetic term gives 

$$-\frac{\psi_{j+1}-2\psi_j+\psi_{j-1}}{\Delta x^2}=\frac{2\psi_j}{\Delta x^2}-\frac{\psi_{j+1}}{\Delta x^2}-\frac{\psi_{j-1}}{\Delta x^2}.$$ 

Therefore, the discrete Hamiltonian entries are 

$$H_{j,j}=\frac{2}{\Delta x^2}+V_j,$$ and $$H_{j,j+1}=H_{j,j-1}=-\frac{1}{\Delta x^2}.$$

This is why the Hamiltonian is tridiagonal for Dirichlet boundary conditions. The main diagonal contains the local kinetic contribution plus the potential, and the off-diagonals contain the nearest-neighbour kinetic coupling.

A very important implementation note is that the same Hamiltonian convention must be used consistently in every comparison. The main Crank--Nicolson step uses the convention $H=-\partial_x^2+V$. Some spectral helper sections in the code use a Hamiltonian written like $H=-\frac{1}{2}\partial_x^2+V$. Both conventions are valid in dimensionless units, but they are not identical. Therefore, analytical or spectral comparisons must use the same coefficient in front of the second derivative as the numerical Crank--Nicolson evolution.

---

## 7. Time discretization

Time is discretized in fixed steps, 

$$t^n=n\Delta t,$$ 

where $\Delta t$ is the timestep. The numerical solution at time step $n$ is the vector 

$$\psi^n=(\psi_0^n,\psi_1^n,\ldots,\psi_{N-1}^n).$$

The choice of $\Delta t$ affects accuracy. Crank--Nicolson is stable in a strong sense, but this does not mean that any timestep gives accurate physics. A very large timestep can still produce inaccurate phases, poor agreement with analytical solutions, and wrong scattering behaviour.

Therefore, convergence tests with different $\Delta x$ and $\Delta t$ are necessary. The repository includes convergence plots for dense Dirichlet, dense periodic and Thomas Dirichlet cases.

---

## 8. Explicit FTCS and why it is not the main method

A simple explicit method would update the wave function using the Hamiltonian evaluated only at the old time level. Schematically, one could write $$\psi^{n+1}=\psi^n-i\Delta t H\psi^n.$$ This is easy to implement because it does not require solving a linear system.

However, this type of explicit method is not ideal for the Schrödinger equation. The exact evolution is unitary, meaning that the norm should be conserved. A basic explicit scheme generally does not preserve unitarity and can introduce artificial growth or damping.

This is why the project focuses on Crank--Nicolson. The cost per timestep is higher because a linear system must be solved, but the qualitative behaviour is much better for quantum dynamics.

---

## 9. Crank--Nicolson method

The Schrödinger equation can be written abstractly as $$i\frac{\partial\psi}{\partial t}=H\psi.$$ Equivalently, $$\frac{\partial\psi}{\partial t}=-iH\psi.$$

The Crank--Nicolson method averages the right-hand side between the old and new time levels. This gives 

$$\frac{\psi^{n+1}-\psi^n}{\Delta t}=-\frac{i}{2}H(\psi^{n+1}+\psi^n).$$ 

Rearranging gives 

$$\left(I+\frac{i\Delta t}{2}H\right)\psi^{n+1}=\left(I-\frac{i\Delta t}{2}H\right)\psi^n.$$

This is written in the code and documentation as 

$$A\psi^{n+1}=B\psi^n.$$

The matrices are 

$$A=I+\frac{i\Delta t}{2}H,$$ and $$B=I-\frac{i\Delta t}{2}H.$$

At each timestep, the code first computes the right-hand side `b = B @ psi`, and then solves `A @ psi_new = b`.

The method is second order in time and, for a Hermitian Hamiltonian, it is unitary in exact arithmetic. This is the key reason for using it.

---

## 10. Matrix coefficients used in the code

The code defines $$r=\frac{i\Delta t}{2\Delta x^2}.$$ For Dirichlet boundary conditions, only the interior points are evolved. The endpoints are fixed to zero.

Using the Hamiltonian $H=-\partial_x^2+V$, the matrix coefficients become:

For matrix $A$, the main diagonal is 

$$A_{j,j}=1+2r+\frac{i\Delta t}{2}V_j,$$

and the off-diagonals are $$A_{j,j+1}=A_{j,j-1}=-r.$$

For matrix $B$, the main diagonal is 

$$B_{j,j}=1-2r-\frac{i\Delta t}{2}V_j,$$

and the off-diagonals are $$B_{j,j+1}=B_{j,j-1}=r.$$

This is exactly what appears in the `CN_step` function. For Dirichlet conditions, the arrays `mainA`, `offA`, `mainB`, and `offB` are created from these expressions.

This is useful because it shows that the numerical code is not arbitrary: each coefficient comes directly from the finite-difference Hamiltonian and the Crank--Nicolson formula.

---

## 11. Dirichlet boundary conditions

Dirichlet boundary conditions impose 

$$\psi(x_{\min},t)=0,\qquad \psi(x_{\max},t)=0.$$ 

Physically, this corresponds to hard walls at the ends of the domain. If a wave packet reaches one of these boundaries, it reflects.

In the implementation, the boundary values are not solved as unknowns. Instead, the code extracts the interior part using `psi_i = psi[1:-1]`. The Crank--Nicolson system is built only for the $N-2$ interior points. After solving, the full wave function is reconstructed with zeros at the endpoints.

This is why the Dirichlet matrix has size $(N-2)\times(N-2)`. It is also why the matrix is tridiagonal: each interior point only couples to its immediate neighbours.

Dirichlet boundaries are appropriate for infinite wells and for finite computational boxes where we want hard-wall behaviour. They are also useful when comparing with the Thomas algorithm because the matrix structure is simple.

---

## 12. Periodic boundary conditions

Periodic boundary conditions impose $$\psi(x_{\min},t)=\psi(x_{\max},t).$$ Physically, this means that the spatial domain behaves like a ring. A wave packet leaving the right side reappears from the left side, and vice versa.

In the matrix, this requires coupling the first and last grid points. The code adds corner terms:

```python
A[0, -1] = -r
A[-1, 0] = -r
B[0, -1] = r
B[-1, 0] = r
```

This changes the matrix from tridiagonal to cyclic tridiagonal. It still has a simple structure, but it is not tridiagonal in the strict Thomas-algorithm sense because of the corner couplings.

For this reason, the current code solves the periodic case using dense linear algebra with `np.linalg.solve`. A modified cyclic Thomas algorithm or Sherman--Morrison method could be implemented later, but the basic Thomas algorithm should not be directly applied to the periodic matrix.

Periodic boundaries are useful because they avoid hard-wall reflection and test a different physical problem.

---

## 13. Dense solver

The dense solver uses `np.linalg.solve(A, b)` to solve the full linear system. This is general and easy to use. It works for both Dirichlet and periodic matrices.

The disadvantage is computational cost. A dense matrix stores all entries, including many zeros, and dense solving becomes expensive for large $N$. For a one-dimensional finite-difference Hamiltonian, most entries are zero, so using a dense solver wastes structure.

However, the dense solver is very useful as a reference implementation. It is simple and robust, and it helps verify that the Thomas solver gives the same result in the Dirichlet case.

---

## 14. Thomas algorithm

The Thomas algorithm is a direct solver for tridiagonal systems. A tridiagonal matrix has nonzero entries only on the lower diagonal, main diagonal and upper diagonal.

A tridiagonal system has the form $$a_i x_{i-1}+b_i x_i+c_i x_{i+1}=d_i.$$ The Thomas algorithm solves it in two stages:

1. Forward elimination: remove the lower diagonal terms.
2. Back substitution: solve for the unknowns from the end of the system backwards.

The cost scales as $$O(N),$$ which is much better than a dense solver for large systems.

In the code, the function `thomas(a, b, c, d)` receives the lower diagonal, main diagonal, upper diagonal and right-hand side. The helper `solveThomas(A, b)` extracts the diagonals from the matrix and calls `thomas`.

Thomas is only used for the Dirichlet case because that matrix is truly tridiagonal. It is not directly used for the periodic case because periodic boundary conditions create corner couplings.

---

## 15. Why Crank--Nicolson plus Thomas is a good combination

Crank--Nicolson gives the correct numerical structure for quantum evolution: stable, time-centred and norm-preserving. The price is that we must solve a linear system at every timestep.

For a one-dimensional finite-difference Hamiltonian with Dirichlet boundaries, the matrix is tridiagonal. Therefore, Thomas solves the Crank--Nicolson system very efficiently.

This combination is ideal for this project:

- Crank--Nicolson gives good quantum time evolution.
- The finite-difference Hamiltonian gives a tridiagonal matrix.
- Thomas exploits the tridiagonal structure.
- The dense solver provides a comparison/reference.

The physics should be the same for dense and Thomas solvers in the Dirichlet case. The difference is computational efficiency.

---

## 16. Norm conservation diagnostic

The norm is computed numerically as $$N(t)=\int |\psi(x,t)|^2dx\approx\sum_j |\psi_j(t)|^2\Delta x.$$ In the code, the norm check uses `np.sum(np.abs(psi)**2) * dx`.

For exact Schrödinger evolution, $N(t)$ should be constant. For a normalized initial state, it should remain close to 1.

The function `comprobar_norma` evolves the system for several valid combinations:

- Dirichlet + dense solver,
- Dirichlet + Thomas solver,
- periodic + dense solver.

The periodic + Thomas case is not physically included as a valid simple Thomas case because the matrix is not strictly tridiagonal.

The plot of $N(t)$ is one of the most important validation results. If Crank--Nicolson is implemented correctly, the norm should remain nearly constant, with only small deviations caused by finite precision, boundary effects or optional renormalization choices.

A strong norm drift would indicate a serious numerical problem.

---

## 17. Renormalization in some sections

Some later parts of the code renormalize the wave function after each step using `psi /= sqrt(trapz(abs(psi)**2,x))`. This is useful to prevent small accumulated numerical errors from affecting plots or reflection/transmission estimates.

However, renormalization should be interpreted carefully. If the method is correctly unitary, it should not need strong renormalization. Therefore, norm conservation should first be tested without hiding errors. After validation, light renormalization can be acceptable for long visual runs or comparison plots.

In a report, it is good to distinguish between:

- norm conservation as a diagnostic, where we check whether the method conserves norm naturally;
- renormalization as a practical plotting or long-run stabilization step.

---

## 18. Analytical reference: infinite well

The infinite well is a useful validation case because the analytical eigenstates are known. If the well extends over a length $L$, the eigenfunctions have sinusoidal shape and the energies are proportional to $n^2$.

In the code, the function `solucion_analitica_pozo` constructs a state of the form 

$$\psi_n(x,t)=\phi_n(x)e^{-iE_nt},$$ 

where 

$$\phi_n(x)=\sqrt{\frac{2}{L}}\sin\left(\frac{n\pi x}{L}\right)$$ 

after shifting the coordinate into the well interval.

If the initial state is one eigenstate, the probability density should not change shape in time because the only time dependence is a complex phase. Therefore, comparing Crank--Nicolson against this analytical solution is a strong test of phase accuracy and boundary handling.

The code stores both numerical and analytical frames and computes an absolute error such as `np.max(np.abs(psi_num - psi_exac))`.

---

## 19. Analytical behaviour: free Gaussian spreading

For a free particle, a Gaussian packet does not simply translate rigidly. It spreads because it is a superposition of different wave numbers. Each wave-number component has a different energy and accumulates phase at a different rate.

The code includes a function `psi_libre_analitica` for the analytical free-packet evolution. This provides a reference for the numerical Crank--Nicolson solution.

The packet width is also tracked. The expected qualitative behaviour is that the width increases with time. This confirms that the simulation is resolving the dispersive nature of the Schrödinger equation.

In the repository, this is connected to the figure `gaussian_packet_spreading_width.png`.

---

## 20. Spectral reference for wells and barriers

For a general finite well or barrier, a simple closed-form time-dependent solution is not always convenient. The code therefore constructs a spectral reference by diagonalizing a discrete Hamiltonian.

The idea is that if the Hamiltonian eigenvectors are $\phi_n$ with eigenvalues $E_n$, an arbitrary initial state can be expanded as 

$$\psi(x,0)=\sum_n c_n\phi_n(x),$$ 

where 

$$c_n=\langle\phi_n|\psi(0)\rangle.$$ 

The time evolution is then 

$$\psi(x,t)=\sum_n c_n\phi_n(x)e^{-iE_nt}.$$

In the code, this is implemented by `np.linalg.eigh(H)`, then projecting the initial state onto each eigenvector, and finally summing the time-evolved modes.

This is called an analytical comparison in the plotting labels, but more precisely it is a spectral numerical reference. It is independent from the Crank--Nicolson time-stepping, so it is still useful for validation.

Again, the Hamiltonian convention must be consistent. If the spectral Hamiltonian uses $-\frac{1}{2}\partial_x^2+V$ while Crank--Nicolson uses $-\partial_x^2+V$, the comparison will contain a scaling mismatch.

---

## 21. Reflection and transmission

For scattering problems, the initial packet is placed to one side of a localized potential and moves toward it. After interaction, part of the wave packet may be reflected and part may be transmitted.

The reflected probability is computed by integrating probability density on the left side of the scattering region, 

$$R(t)=\int_{x<x_L}|\psi(x,t)|^2dx.$$ 

The transmitted probability is computed on the right side, 

$$T(t)=\int_{x>x_R}|\psi(x,t)|^2dx.$$

In the code, this is done using masks such as `x < x_left` and `x > x_right`, then applying `np.trapz` to $|\psi|^2$ in those regions.

After the packet has fully interacted with the potential and moved away, one expects approximately $$R+T\approx1,$$ assuming the integration regions cover all outgoing probability and the boundaries have not introduced additional effects.

This diagnostic is physically important because it connects the numerical wave function with measurable scattering probabilities.

---

## 22. Boundary effects in scattering

When studying reflection and transmission, the computational domain must be large enough. The wave packet should start far enough from the potential, and the boundaries should not affect the scattering event too early.

If the domain is too small, reflected waves from the boundary can mix with the physical reflected/transmitted waves from the potential. Then $R(t)$ and $T(t)$ become difficult to interpret.

Dirichlet boundaries are hard walls, so they can create artificial reflections if the packet reaches them. Periodic boundaries can wrap the packet around the domain, which is also not appropriate for all scattering interpretations.

Therefore, scattering simulations require careful choices of:

- domain length $L$,
- initial packet position $x_0$,
- packet width $\sigma$,
- wave number $k_0$,
- potential width and height,
- total simulation time,
- regions used to calculate $R$ and $T$.

---

## 23. Convergence studies

The repository includes convergence diagnostics through the function `estudio_convergencia_CN`. The idea is to run the same problem with different values of $\Delta x$ and $\Delta t$, and compare the numerical solution to a reference.

A good method should improve as the grid is refined. For Crank--Nicolson with centred finite differences, the expected formal accuracy is second order in time and second order in space under suitable conditions.

However, the observed error can depend on many things:

- the reference solution used,
- the boundary conditions,
- the total simulation time,
- whether the packet reaches the boundary,
- the smoothness of the potential,
- the compatibility of the initial state with the grid,
- the Hamiltonian convention.

The convergence plots are therefore important because they show whether the code behaves consistently when resolution is changed.

---

## 24. Computational-cost comparison

The code also measures computational time for different solvers and grid sizes. This is important because the dense solver and Thomas solver can give the same physical answer but with different costs.

The dense solver is easier to implement and works for more general matrices, but it does not exploit sparsity. The Thomas algorithm is specialized, but for tridiagonal systems it is much faster and scales better with $N$.

The expected behaviour is:

- for small $N$, the difference may not be very important;
- for large $N$, Thomas should become much more efficient for Dirichlet problems;
- periodic dense solving is more expensive because the current implementation uses a full matrix solve.

This demonstrates a key numerical lesson: the mathematical structure of the discretized problem should guide the solver choice.

---

## 25. Main implementation workflow

The main numerical workflow can be summarized as:

```python
# 1. Choose spatial grid
L = 10
N = 800
x = np.linspace(-L, L, N)
dx = x[1] - x[0]

# 2. Choose timestep and number of steps
dt = 0.1
steps = 2000

# 3. Define potential
V = V_func(x)

# 4. Build initial Gaussian packet
psi0 = np.exp(-(x-x0)**2/(2*sigma**2)) * np.exp(1j*k0*x)
psi0 /= np.sqrt(np.trapz(np.abs(psi0)**2, x))

# 5. Evolve using Crank--Nicolson
psi = psi0.copy()
for n in range(steps):
    psi = CN_step(psi, V, dx, dt, resolver="thomas", cc_tipo="dirichlet")

# 6. Analyse probability density
rho = np.abs(psi)**2
norm = np.sum(rho)*dx
```

This is the essential structure of the project.

---

## 26. Crank--Nicolson step in algorithm form

The function `CN_step` does the following:

```python
def CN_step(psi, V, dx, dt, resolver="normal", cc_tipo="dirichlet"):
    # Compute r = i dt / (2 dx^2)
    r = 1j * dt / (2 * dx * dx)

    # If Dirichlet:
    #   Use only interior points psi[1:-1]
    #   Build tridiagonal A and B
    #   Compute b = B @ psi_interior
    #   Solve A psi_new = b using dense or Thomas
    #   Reconstruct full psi with endpoints equal to zero

    # If periodic:
    #   Use all points
    #   Build A and B with corner couplings
    #   Compute b = B @ psi
    #   Solve using dense solver

    return psi_new
```

The important idea is that the same mathematical method is used in both boundary conditions, but the matrix structure changes.

---

## 27. Thomas algorithm in algorithm form

The Thomas solver works as:

```python
def thomas(a, b, c, d):
    # a: lower diagonal
    # b: main diagonal
    # c: upper diagonal
    # d: right-hand side

    # Forward elimination
    for i in range(1, n):
        m = a[i-1] / b[i-1]
        b[i] = b[i] - m*c[i-1]
        d[i] = d[i] - m*d[i-1]

    # Back substitution
    x[-1] = d[-1] / b[-1]
    for i in reversed(range(n-1)):
        x[i] = (d[i] - c[i]*x[i+1]) / b[i]

    return x
```

This is why Thomas is efficient: it never constructs or manipulates a full dense matrix during the elimination.

---

## 28. Output figures and diagnostics

The repository is organized with figures showing different aspects of the simulation:

- `figures/free_packet/`: free Gaussian propagation and spreading.
- `figures/boundary_conditions/`: comparison between boundary conditions and solvers.
- `figures/validation/`: norm conservation and convergence checks.
- `figures/analytical_comparison/`: comparisons with analytical or spectral references.
- `figures/scattering/`: reflection and transmission for wells and barriers.
- `figures/performance/`: computational-cost comparisons.

This organization is useful because each figure group corresponds to a different validation question:

- Does the packet move correctly?
- Does the norm stay constant?
- Do dense and Thomas solvers agree?
- Does the error decrease with resolution?
- Is the Thomas algorithm faster?
- Do scattering probabilities behave physically?

---

## 29. Interpretation of free-packet results

For the free particle, the expected behaviour is that the wave packet propagates in the direction set by $k_0$ and spreads over time. The spreading is not a numerical error; it is a physical consequence of the wave packet being a superposition of momenta.

Under Dirichlet boundary conditions, if the packet reaches the edge of the box, it reflects because the boundary acts as a hard wall.

Under periodic boundary conditions, if the packet reaches the edge, it reappears from the other side. This is not the same physical problem as Dirichlet boundaries, so the two animations should not be interpreted as two solvers of the same physical setup. They represent different boundary physics.

---

## 30. Interpretation of solver comparison

For the Dirichlet case, dense solving and Thomas solving should produce the same evolution. The Crank--Nicolson matrices are the same; only the method of solving the linear system changes.

Therefore, if the dense and Thomas animations overlap, this validates the Thomas implementation.

If they disagree significantly, possible causes include:

- incorrect extraction of diagonals,
- wrong boundary treatment,
- applying Thomas to a non-tridiagonal matrix,
- inconsistent normalization,
- programming error in matrix construction.

In the current structure, Thomas is appropriate for Dirichlet but not for the simple periodic implementation.

---

## 31. Interpretation of norm-conservation results

The norm-conservation plot checks whether the numerical method respects probability conservation. For a good Crank--Nicolson implementation, the norm should remain very close to the initial value.

Small numerical deviations are acceptable because of floating-point roundoff and finite-grid effects. However, large monotonic drift would be a serious problem.

This diagnostic is especially important because a visually plausible animation can still be wrong if the norm is not conserved. Therefore, norm conservation is one of the strongest tests in the project.

---

## 32. Interpretation of infinite-well comparison

In an infinite well, analytical eigenstates are known. If the initial wave function is an eigenstate, the probability density should remain stationary. Only the complex phase changes.

The comparison between Crank--Nicolson and the analytical solution checks that the numerical method correctly reproduces this phase evolution and does not distort the density.

The absolute error plot gives a quantitative measure. A smaller error indicates better agreement. The error should decrease when $\Delta x$ and $\Delta t$ are refined, assuming the reference solution uses the same physical convention.

---

## 33. Interpretation of harmonic-oscillator comparison

For the harmonic oscillator, a localized Gaussian packet should oscillate in the confining potential. This is consistent with Ehrenfest's theorem: the expectation value of position follows the classical oscillator motion for suitable coherent-state-like packets.

The harmonic oscillator test is useful because it checks bound motion rather than free propagation or hard-wall reflection.

If the numerical packet centre follows the analytical centre reasonably well, the simulation is capturing the expected dynamics. Differences can come from finite grid size, packet not being an exact coherent state, boundary effects, timestep error or inconsistent Hamiltonian scaling.

---

## 34. Interpretation of scattering results

For the finite well and rectangular barrier, the wave packet splits into reflected and transmitted parts. The functions `calcular_reflexion_transmision`, `graficar_instantes_pozo`, and `graficar_instantes_barrera` compute and plot this behaviour.

For a barrier, increasing the barrier height should generally increase reflection and reduce transmission, although quantum tunnelling can still allow nonzero transmission.

For a well, reflection can also occur because the wave packet sees a sudden change in potential. The well does not simply absorb the packet; it creates interference and partial reflection/transmission.

A physically meaningful scattering result should satisfy approximately $R+T\approx1$ after the packet has separated, assuming the left and right integration regions include all probability and boundaries have not interfered.

---

## 35. Limitations and important cautions

The project is a strong educational implementation, but some limitations should be stated clearly.

First, the equation is written in dimensionless units. Therefore, the numerical values of $x$, $t$, $V$, $k_0$ and $\sigma$ are not automatically physical SI quantities.

Second, boundary conditions define the physical problem. Dirichlet and periodic simulations are not expected to give the same behaviour near the edges.

Third, Thomas should only be used for strictly tridiagonal systems. The periodic matrix has corner couplings, so a simple Thomas algorithm is not directly valid there.

Fourth, analytical and spectral comparisons must use the same Hamiltonian convention as the Crank--Nicolson evolution. In particular, the factor in front of the kinetic term must be consistent.

Fifth, renormalizing after every step can hide norm-conservation errors. It is useful for stable plots, but the method should first be checked without relying on renormalization.

Sixth, dense solvers are convenient but become expensive for large grids. A sparse or cyclic-tridiagonal solver would be a natural improvement for larger simulations.

---

## 36. Final summary

This project solves the one-dimensional time-dependent Schrödinger equation for Gaussian wave packets under different potentials and boundary conditions. The physical quantity of interest is the probability density $|\psi(x,t)|^2$, and the most important conservation law is conservation of total probability.

The wave function is discretized on a uniform grid. The second derivative is approximated using a centred finite-difference stencil, which produces a nearest-neighbour Hamiltonian matrix. Time evolution is performed with the Crank--Nicolson method, written as $A\psi^{n+1}=B\psi^n$.

Crank--Nicolson is chosen because it is stable, second-order accurate and well suited to unitary quantum evolution. For Dirichlet boundary conditions, the matrices are tridiagonal, so the Thomas algorithm can solve the linear system efficiently. For periodic boundary conditions, the first and last grid points are coupled, so the current implementation uses dense solving.

The code studies free propagation, infinite wells, harmonic confinement, finite wells and rectangular barriers. The diagnostics include norm conservation, convergence/error plots, dense-versus-Thomas solver comparisons, computational-cost measurements, analytical or spectral reference comparisons, Gaussian packet spreading, and reflection/transmission probabilities.

Overall, the repository demonstrates a complete numerical workflow for one-dimensional computational quantum mechanics: define the physical problem, discretize it, evolve with a norm-preserving method, compare solvers, validate against references, and interpret the resulting wave-packet dynamics.


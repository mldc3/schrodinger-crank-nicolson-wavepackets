# Sources and Notes

## Code provenance

The main implementation is preserved in:

- `src/schrodinger_crank_nicolson.py`

This repository is a cleaned portfolio version of original coursework code and figures. The numerical logic and original physical parameter choices were intentionally preserved.

## Coursework and personal-notes foundation

The project content is based on coursework and personal notes in computational physics / modelling, especially:

- finite-difference spatial discretization,
- PDE time evolution workflows,
- explicit versus implicit schemes,
- FTCS limitations for Schrödinger dynamics,
- Crank–Nicolson time integration,
- tridiagonal linear systems,
- Thomas algorithm implementation,
- 1D time-dependent Schrödinger equation modeling,
- Gaussian wave-packet initialization and propagation,
- norm-conservation diagnostics,
- analytical/spectral validation,
- scattering diagnostics via reflection/transmission estimates.

## Original report reference

If present in the repository, the original coursework report is stored at:

- `notes/original_course_report/practica7parteII.pdf`

## Organization note

The current layout separates:

- source code (`src/`),
- technical documentation (`docs/`),
- generated figures (`figures/...`), and
- archival/report material (`notes/...`).

## Scope note

This polish and reorganization does not claim new simulations, new physics, or modified numerical algorithms; it improves presentation, structure, and reproducibility of existing work.

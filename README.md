# Bloch-Redfield Educational Repository

This mini-project provides pedagogical Jupyter notebooks illustrating the Bloch-Redfield (BR) master equation using QuTiP's `brmesolve` and comparisons to Lindblad (`mesolve`) where appropriate.

## Goals
- Show simple canonical open quantum system examples.
- Visualize dynamical maps and identify when BR breaks CP (negative density matrix eigenvalues).
- Demonstrate the secular approximation: Redfield vs Secular Redfield vs Lindblad.
- Explore multiple bath spectral densities (`DrudeLorentzEnvironment`, `OhmicEnvironment`, `UnderDampedEnvironment`).

## Notebooks
1. `01_qubit_dephasing.ipynb`: Pure dephasing qubit: analytical vs BR vs Lindblad.
2. `02_qubit_relaxation_drude_lorentz.ipynb`: Energy relaxation & temperature dependence.
3. `03_harmonic_oscillator_ohmic.ipynb`: Damped quantum oscillator & approach to thermal state.
4. `04_coupled_qubits_secular_breakdown.ipynb`: Breakdown of secular approximation when splittings comparable to rates.
5. `05_negative_eigenvalues_explorer.ipynb`: Scan parameter regimes to locate positivity violations.

## Structure
`src/` contains small reusable helpers (bath factory, diagnostics, plotting helpers).
`tests/` minimal sanity tests (import & a short evolution) to ensure environment reproducibility.

## Install
Create environment (example):
```
conda create -n redfield-ed python=3.11 -y
conda activate redfield-ed
pip install qutip numpy scipy matplotlib
```
(Optional) Install this package editable:
```
pip install -e .
```

## License
MIT (add if desired).

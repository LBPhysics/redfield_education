# Bloch-Redfield Educational Repository

This mini-project provides pedagogical Jupyter notebooks illustrating the Bloch-Redfield (redfield) master equation using QuTiP's `brmesolve` and comparisons to Lindblad (`mesolve`) where appropriate.

## Goals
- Show simple canonical open quantum system examples.
- Visualize dynamical maps and identify when redfield breaks CP (negative density matrix eigenvalues).
- Demonstrate the secular approximation: Redfield vs Secular Redfield vs Lindblad.
- Explore multiple bath spectral densities (`DrudeLorentzEnvironment`, `OhmicEnvironment`, `UnderDampedEnvironment`).

## Notebooks
1. `01_qubit_dephasing.ipynb`: Pure dephasing qubit: analytical vs redfield vs Lindblad.
2. `02_qubit_relaxation_drude_lorentz.ipynb`: Energy relaxation & temperature dependence.
3. `03_harmonic_oscillator_ohmic.ipynb`: Damped quantum oscillator & approach to thermal state.
4. `04_coupled_qubits_secular_breakdown.ipynb`: Breakdown of secular approximation when splittings comparable to rates.
5. `05_negative_eigenvalues_explorer.ipynb`: Scan parameter regimes to locate positivity violations.

## Structure
`src/` contains small reusable helpers (bath factory, diagnostics, plotting helpers).
`tests/` minimal sanity tests (import & a short evolution) to ensure environment reproducibility.

## Install

### Option 1: Using conda environment file (Recommended)
```bash
# Create and activate environment from environment.yml
conda env create -f environment.yml
conda activate redfield-education

# Install this package in editable mode (optional)
pip install -e .
```

### Option 2: Manual environment creation
```bash
# Create environment manually
conda create -n redfield-ed python=3.11 -y
conda activate redfield-ed
pip install qutip numpy scipy matplotlib jupyter pytest black ruff tqdm

# Install this package in editable mode (optional)
pip install -e .
```

### Updating the environment
```bash
# Update environment from environment.yml
conda env update -f environment.yml --prune
```

## License
MIT (add if desired).

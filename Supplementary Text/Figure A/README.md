# Supplementary Figure A: Latin hypercube sensitivity

This directory contains Latin hypercube sampling runs exploring sensitivity to release parameters.

## Contents
- `main_Latincube_YLE.R` – Runs Latin hypercube simulations for YLE scenarios.
- `generate_YLE_inheritance_cube.R` – Builds the inheritance cube used in the runs.
- `FigA_Sensitivity_BRT.R` – Processes outputs and generates the supplementary figure.
- `mgdriveYLE_sweep_Latin4/` – Stored simulation outputs from the Latin hypercube sweep.

Run the main script from this directory to regenerate simulations, then use `FigA_Sensitivity_BRT.R` to produce the figure.

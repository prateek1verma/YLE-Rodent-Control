# Figure 1 data generation

This directory runs the MGDrivE sweeps that underpin Figure 1.

## Key scripts
- `main_YLE_vary_rel_fy_target_fertility_haploinsufficient.R` – Configures and executes YLE simulations targeting female fertility gene across different release rates while keeping fitness cost to zero.
- `generate_YLE_inheritance_cube.R` – Builds the inheritance cube used by the simulations.

## Outputs
Simulation outputs are organized under `mgdriveYLE_sweep_10yr_release_target_fertility/` with nested folders for each parameter set.

## Usage
Run the main script in R from this directory. The script loads MGDrivE (MouseGD), sets parameter grids, and saves results into the structured output folders. Change the parameters in the main script to run simulations for the desired scenario.

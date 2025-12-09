# Figure 1: Fertility-targeting YLE sweeps

This folder contains the simulations and plotting scripts that generate Figure 1.

## Structure
- `Data_generator_Fig1/` – Runs MGDrivE (MouseGD) sweeps for fertility-targeting YLE releases and stores simulation outputs in csv files.
- `Fig_generator_Fig1/` – Creates the individual panels (e.g., `Fig1b`, `Fig1c`) from processed data.
- `Figure_panel/` – Final assembled figure assets.

## Reproducing the simulations
1. Ensure MGDrivE (MouseGD) is installed (see root README).
2. From `Data_generator_Fig1/`, run `main_YLE_vary_rel_fy_target_fertility_haploinsufficient.R` to regenerate sweep outputs. Change parameters as required, for eg. number of runs to 100 runs. The `generate_YLE_inheritance_cube.R` helper constructs the inheritance cube for each run.
3. Generated data are stored in the `mgdriveYLE_sweep_10yr_release_target_fertility/` subfolders. The provided output can be reused to skip long simulations.

## Regenerating plots
From `Fig_generator_Fig1/`, run `Fig1b_plot.R` and `Fig1c_plot.R` to recreate the panels saved as PNG files.

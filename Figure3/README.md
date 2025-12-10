# Figure 3: Haplosufficient vs. haploinsufficient targets

This folder houses simulations and plots comparing target effect assumptions for YLE strategies.

## Structure
- `Data_generator_Fig3/` – Runs simulations for both haplosufficient and partially haploinsufficient targets.
- `Fig_generator_Fig3/` – Generates figure panels from summary data.
- `Panel_Figure3/` – Final combined figure assets.

## Reproducing the simulations
Run the scripts in `Data_generator_Fig3/` to regenerate the parameter sweeps. It uses dependent function defined in `generate_YLE_inheritance_cube.R` to create the inheritence cube in the main scripts.

## Regenerating plots
Execute `Fig_generator_Fig3/Fig3_plots.R` to rebuild the consolidated figure using the processed outputs.

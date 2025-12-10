# Figure 5: Multi-patch comparisons

Contains simulations and plotting scripts for comparing YLE, homing, and X-shredder strategies across multiple patches.

## Structure
- `Data_generator_Fig5/` – Runs multi-patch simulations for each strategy.
- `Fig_generator_Fig5/` – Produces the individual panels.
- `Panel_Figure5/` – Final combined figure assets.

## Reproducing the simulations
Use the scripts in `Data_generator_Fig5/` to run the multi-patch parameter sweeps. Helper scripts such as `cube_XShredderY.R`, `generate_YLE_inheritance_cube.R`, and `generate_Homing_inheritance_cube.R` configure the inheritance structures. Main scripts (`main_YLE_multi_patches.R`, `main_Homing_many_patch_*`, `main_Xshredder_many_patch_*`) run the simulations.

## Regenerating plots
Run the plotting scripts in `Fig_generator_Fig5/` (`Fig5a.R`, `Fig5b.R`, `Fig5c.R`, `Fig5d.R`) to recreate the panels before assembling them in `Panel_Figure5/`.

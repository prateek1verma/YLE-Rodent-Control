# Figure 2: Sensitivity analyses

Materials for Figure 2 exploring YLE sensitivity to release and fitness parameters.

## Structure
- `Data_generator_Fig2/` – Runs parameter sweeps and summarises outputs.
- `Fig_generator_Fig2/` – Produces the panel plot from processed data.
- `Figure2_panel/` – Final assembled figure asset.

## Reproducing the simulations
1. Install required R packages (see root README for MGDrivE installation).
2. From `Data_generator_Fig2/`, run `main_YLE_vary_rel_fy_target_fertility_haploinsufficient.R` to execute sweeps. Use `generate_YLE_inheritance_cube.R` for cube setup and `Process_large_data_create_summary_csv_files.R` to summarise results for plotting.

## Regenerating plots
Run `Fig_generator_Fig2/Fig2_plot.R` to rebuild the figure using the summary CSV output generated using data generated with `main_YLE_vary_rel_fy_target_fertility_haploinsufficient.R`.

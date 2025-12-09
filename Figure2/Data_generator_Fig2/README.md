# Figure 2 data generation

Scripts for running the YLE sensitivity sweeps that feed Figure 2.

## Key scripts
- `main_YLE_vary_rel_fy_target_fertility_haploinsufficient.R` – Executes the parameter grid for releases and fitness costs (as reduction in life-span).
- `generate_YLE_inheritance_cube.R` – Builds the inheritance cube required by MGDrivE.
- `Process_large_data_create_summary_csv_files.R` – Aggregates simulation outputs into summary CSV files for plotting.

Run the main script from this directory. Output data are saved locally and can be re-used by the plotting script.

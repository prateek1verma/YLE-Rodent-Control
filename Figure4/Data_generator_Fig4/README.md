# Figure 4 data generation

Runs alternative drive architecture simulations for Figure 4.

## Key scripts
- `main_fsRRDL_para_sweep.R` – Parameter sweeps for fsRDDL.
- `main_SMR_para_sweep.r` – Parameter sweeps for self-limiting releases.
- `generate_YLE_inheritance_cube.R`, `generate_fRIDL_inheritance_cube.R`, `generate_SMR_inheritance_cube.R` – Build inheritance cubes for each drive type.
- `Process_large_data_create_summary_csv_files.R` – Summarises outputs for plotting.

Run the desired sweep scripts in R; outputs are saved locally for downstream plotting.
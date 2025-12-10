# Figure 4: Alternative drive architectures

Contains simulations and plotting scripts for Figure 4 comparing Sterile Male Releases (SMR), fsRDDL, and YLE strategies.

## Structure
- `Data_generator_Fig4/` – Runs parameter sweeps for alternative drive architectures and processes outputs.
- `Fig_generator_Fig4/` – Generates individual panels.
- `Panel_Figure4/` – Final combined figure assets.

## Reproducing the simulations
Use the scripts in `Data_generator_Fig4/` to run the sweeps. The inheritance cube helpers (`generate_YLE_inheritance_cube.R`, `generate_fRIDL_inheritance_cube.R`, `generate_SMR_inheritance_cube.R`) support the different genetic control strategies. `main_fsRRDL_para_sweep.R` and `main_SMR_para_sweep.r` coordinate the parameter grids, and `Process_large_data_create_summary_csv_files.R` summarises outputs.

## Regenerating plots
Run `Fig_generator_Fig4/Fig4d.R` and `Fig_generator_Fig4/Fig4e.R` to build the respective panels in `Panel_Figure4/`.

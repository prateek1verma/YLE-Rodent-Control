# Figure 6 data generation

Simulations supporting Figure 6's multi-patch comparisons and heatmap.

## Key scripts
- `main_contritutive_YLE_multi_patches.R` – Boosted YLE simulations across connected patches. Change the library(MGDrivEmouse2) to whatever is the name of the MouseGD package you have locally installed on your machine.
- `main_contritutive_YLE_vary_rel_xshred_heatmap.R` – Plots the heatmap for Boosted YLE simulations while varying release size and efficiency of X-Shredder. 

Change the library(MGDrivEmouse2) to whatever is the name of the MouseGD package you have locally installed on your machine.

- `generate_boosted_YLE_constitutive_inheritance_cube.R` – Build inheritance structures for boosted YLE.

Run the desired main scripts in R from this directory. Outputs are saved locally for plotting.

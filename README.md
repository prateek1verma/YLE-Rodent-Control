# Y-linked editors for invasive rodent control: A mathematical modeling study

**Authors:** Prateek Verma, Omar S. Akbari, John M. Marshall  
**Status:** In Review / Pre-print  

## 📌 Overview

This repository contains the source code, simulation scripts, and data analysis workflows for the research paper *"Y-linked editors for invasive rodent control: A mathematical modeling study"*.

In this study, we evaluate the potential of **Y-linked genome editors (YLEs)**—a self-limiting genetic biocontrol strategy—to suppress or eliminate invasive mouse populations. We benchmark YLEs against other self-limiting genetic biocontrol tools (such as sterile male releases) and self-sustaining tools (such as Homing Gene Drives and Y-linked X-shredders) using the **MouseGD** stochastic simulation framework.

## 📂 Repository Structure

The repository is organized by Figure to facilitate the reproduction of specific results found in the manuscript:

* **`Figure1/`**: Parameter sweeps and plots for fertility-targeting YLE scenarios.
* **`Figure2/`**: Sensitivity analyses varying key release and fitness parameters.
* **`Figure3/`**: Comparisons of haplosufficient vs. partially haploinsufficient target effects.
* **`Figure4/`**: Alternative drive architectures including fRIDL and fsRDDL.
* **`Figure5/`**: Multi-patch releases contrasting YLE containment with Homing and X-shredder invasion.
* **`Figure6/`**: Analysis of breeding facility requirements and graphical abstract assets.
* **`Supplementary Text/`**: Extended sensitivity analyses and Boosted Regression Tree (BRT) models (Figures A–H).

## 🛠️ Installation & Prerequisites

The analysis is performed in **R**. The simulations rely on the `MGDrivE` framework.

### 1. Install MGDrivE (MouseGD)

This project utilizes the **MouseGD** framework, implemented within the `MGDrivE` package as described by Brown *et al.* (2022). To reproduce the simulations, please install the package directly from the source repository:

**Reference:**

> Brown EA, Eikenbary SR, Landis WG. Bayesian network‐based risk assessment of synthetic biology: Simulating CRISPR‐Cas9 gene drive dynamics in invasive rodent management. *Risk Analysis*. 2022 Dec;42(12):2835-46.

**Installation Command:**

```r
library(devtools)
# Install MGDrivE from the specific GitHub repository
install_github("eabrown2378/MGDrivE/MGDrivE")
```

## License
This project is released under the GPL-3.0 License and may be modified or extended as needed. See the LICENSE terms for details.

## 📧 Contact

For questions regarding the code or model parameters, please contact:

  * **Prateek Verma** - prateekverma@berkeley.edu
  * **John M. Marshall** - john.marshall@berkeley.edu

<!-- end list -->

# Synthetic biology design principles enable efficient bioproduction of Heparosan with low polydispersion index for the biomedical industry
This is the repository of the code we used to generate the results in "Synthetic biology design principles enable efficient bioproduction of Heparosan with low polydispersion index for the biomedical industry" manuscript.

## Matlab Version

We used Matlab version 23.2.0.2485118 (R2023b) Update 6, we did not test if this code works with other Matlab versions.

## Repository structure

Included within this directory are all the scripts required to execute the simulations and optimizations, yielding the results showcased in the manuscript. This repository encapsulates everything essential for reproducing the findings detailed in the paper.

### /optimization

**ODE Model**

The file `model_hep6p.m` contains the ODE model definition. The parameters of the model are implemented in the function `parameters.m`. 

**Files to use MEIGO optimization algorithm**

* File `CostFunction_HEP.m` contains the implementation of the cost function for our multiobjective optimization problem.
* The file `optimize_with_MEIGO.m` is used to call the optimization algorithm used [MEIGO](https://github.com/gingproc-IIM-CSIC/MEIGO64). Check the `EXAMPLES` folder within the MEIGO repository for more info on how to use it. In this case, the file `optimize_with_MEIGO.m` includes all the parameters necessary for the execution with our cost function `CostFunction_HEP.m`.

### /robustness 

Scripts for the robustness analysis.
* `Simulation_and_plot_perturbation.m` script is used to get the different plots for the manuscript to illustrate the results.
* `Optimal_Values_Calculation.mlx` script calculates the uncertainties of the optimum values using the gaussian process.
* `distributuions_optim.mlx` and `distributuions_random.mlx` scripts calculate the optimum and random distributions using values of PDI and Mw for the two sets of data: optimal values from the optimization and random values.
* `level_diagram_plot.m` plot the PDI data against MW, and the precursors concentraitons for Figure 4 of the manuscript.
* `results_cL1_10_cL2_10_cU1_20_cU2_20_date_20240331_154057.mat` has the optimization results in Matlab format. 
* `tablaOutput.mat` compiles the results in a table format.                                              
### /fba

We include all the code used to develop the flux balance analysis in different files:

* The original metabolic SBML model **iDK1463**, written as a MATLAB object `nissle_core.mat`.
* The new metabolic SBML model including the heparosan biosynthesis pathway (metabolites, enzymes, and reactions), written as matlab object `nissle_core_heparosan.mat`.
* The COBRA file `FBA_heparosan.mlx` (Matlab LiveScript) to add the heparosan biochemical reactions, and to perform the FBA.





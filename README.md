# Synthetic biology design principles enable efficient bioproduction of Heparosan with low polydispersion index for the biomedical industry
This is the repository of the code we used to generate the results in "Synthetic biology design principles enable efficient bioproduction of Heparosan with low polydispersion index for the biomedical industry" manuscript.

## Matlab Version

We used Matlab version 23.2.0.2485118 (R2023b) Update 6, we did not test if this code works with other Matlab versions.

## Repository structure

**scripts/** 
 
Included within this directory are all the scripts required to execute the simulations and optimizations, yielding the results showcased in the manuscript. This repository encapsulates everything essential for reproducing the findings detailed in the paper.

**Files to use MEIGO optimization algorithm**

* The file `optmize_with_MEIGO.m` is used to call the optimization algorithm used [MEIGO](https://github.com/gingproc-IIM-CSIC/MEIGO64). Check the `EXAMPLES` folder within the repository for more info on how to use it. In this case, the file `optimize_with_MEIGO.m` includes all the parameters necessary for the execution with our cost function `CostFunction_HEP.m`.
* File `CostFunction_HEP.m` contains the implementation of the cost function for our multiobjective optimization problem.

**ODE Model**

The file `model_hep6p.m` contains the ODE model definition. The file `Main_anthitetic_marker.m` contains the script to run the model, and to obtain simulation results like in the next figure. The parameters of the model are implemented in the function `parameters.m`. 

**scripts/for_plot** 

* Scripts necessary to obtain the plots of the manuscript.




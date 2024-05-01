% Optimization of PDI and Mw of Hepraosan using ODE model for the genetic
% and metabolic circuits and GP to model the relationship between the PDI
% Mw and the precursors.
% Parameters: structure contains all rates and constants
% Updated 28/04/2024 by Yadira Boada, Alejandro Vignoni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear mex;
clear all;
close all;
dbstop if error

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='CostFunction_HEP'; %mfile containing the objective function
 
problem.x_L = [0.1 0.1 0.1 ]; %Lower bound for parameters
problem.x_U = [ 10 10 4.8 ];  %Upper bound for parameters
problem.x0 = [7.1273    1.7276    1.6621]; %Initial guess parameters

problem.neq = 0; % Number of equalities
problem.c_L = [10 10  0 0 0 0 -7.5 -7.5]; %Lower bound constraint
problem.c_U = [20 20  7.5 7.5 7.5 7.5  0  0];  %Upper bound constraint

opts.maxeval=500;
opts.ndiverse=25;
opts.local.solver='dhc';
opts.local.finish='dhc';
opts.local.iterprint=1;
opts.inter_save = 1;
opts.plot = 1; %Plot convergence at the end

%========================= END OF PROBLEM SPECIFICATIONS =====================

Results=MEIGO(problem,opts,'ESS',0);

% Get the current date and time as a string
current_datetime = datestr(now, 'yyyymmdd_HHMMSS');

% Create a filename with the numeric variable and date/time
filename = sprintf('results_cL1_%d_cL2_%d_cU1_%d_cU2_%d_date_%s.mat',problem.c_L(1),problem.c_L(2),problem.c_U(1),problem.c_U(2) , current_datetime);
save(filename, 'Results');



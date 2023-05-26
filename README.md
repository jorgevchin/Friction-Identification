# Friction-Identification
This repository contains the source code and data files for the manuscript entitled:

"A Novel Method for Obtaining the Friction Models of Direct-Drive Robot Manipulators: Characterization and Quantitative Comparison".

The manuscript is currently under review at IEEE Latin America Transactions.

The results were generated using MATLAB R2021a (9.10.0.1602886) 64-bit (win64) February, 23021

## Folder description

- Comparacion
  Simulation.slx: Simulink model that simulates the dd2dof robot dynamic model control.
  Simulations.m: Script that automatically generates the simulation data for comparison.
  Comparison.m: Script that compares the simulation results with experimental results.
  Files with extension .mat contains the simulation results.

- Test 1
  Test1.m: Script that calculates the static friction from the results of test 1.

- Test 2
  Test2.m: Script that calculates the Viscous and Coulomb friction coefficients and the fitting parameters in the Stribeck effect 

- Test 3
  Test3.m: Calculates the stiffness and damping factors in the Dahl and LuGre models
  
The files with extension .t contains the experimental results


## Instructions for reproducing results
1) Perform the experimental tests as described in the manuscript.
2) Run the Test#.m scripts in the Test 1,2,3 folders in order to estimate the friction parameters. These parameters are then introduced into Simulation.slx. 
3) Run Simulations.m to generate all the necessary files for comparison. 
4) Run Comparision.m to obtain a quantitative comparison between the simulations with each friction model and the experimental data.

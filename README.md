# MigrationModeling
Writen by Joanna Bieri and Christine Sample

This code is based on the general network framework as described in the paper:

A general modeling framework for describing spatio-temporal population dynamics

Christine Sample* , John Fryxell, Joanna Bieri, Paula Federico, Julia Earl, Ruscena Wiederholt, Tyler Flockhart, Sam Nicol, Jay E. Diffendorfer, Brady J. Mattsson, Wayne E. Thogmartin, Richard A. Erickson, and D. Ryan Norris.

Version 1.2 of the code as of 5/7/18

Updates include major changes to the structure of the code to accomodate the new more general matrix version of the CR calculation. This update allows the Cr calculation to be easily implemented into the network code, but changes how we think about the f-update function. Rather than output of numbers of individuals, this function now must output node update "growth or decay" rate values. Then the time step update in the network code is a matrix multiplication (descrpibed in Sample et al CR Kr paper) This version of the code is NOT COMPATABLE with 1.1, since overall handling of data: matrix vs. list representation, has been changed!!!



Version 1.1 of the code as of 8/3/17

Requires the following R libraries:
library(XLConnect)
** Updated to remove error on Macs that causes outputs to fail
** Updated tutorial - see tutorial for more information on running the code.



NetworkSetup.R - Reads in the Network Parameters which are stored in .xlsx files. This code depends on the format of the .xlsx files. The formats must be followed carefully in order for this to work.

NetworkSimulation.R - Executes the general network model as described in the paper. Depends on:
*** NetworkSetup.R (reading of parameters in the .xlsx files) 
*** SpeciesFunctions.R (correct definition of f_(i,t), p_(ij,t) and s_(ij,t))
*** Original parameters as set up in the SpeciesSimulation.R file.

NetworkOutputs.R - Gives a simple graphical and numerical output of the simulation results.


* To run the network model for your species,
  - define species-specific parameters in the input file(s):  ./SIMNAME/network_inputs_NETNAME.xlsx
  - define species-specific functions f_(i,t), p_(ij,t) and s_(ij,t) in SpeciesFunctions.R
  - define species-specific data (SIMNAME, NETNAME, seasons, and num_nodes) in SpeciesSimulation.R
  - run SpeciesSimulation.R

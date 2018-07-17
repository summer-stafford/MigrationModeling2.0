##################
## NETWORK CODE ##
##################
library(XLConnect)


#########################################
### SET SPECIES SPECIFIC NETWORK INFO ###
#########################################

### MONARCH EXAMPLE ###

# IMPORTANT NOTE #
# Users should specify network model functions f_(i,t), p_(ij,t), and s_(ij,t) in the file:
# NetworkFunctions.R, even if these are constants they must still be set as such.

### User specified information specific to the simulation ###

SIMNAME <- "Baseline1" # Specifies which baseline folder 
# This is the subfoulder where the network data, outputs, and .RData file should are saved
# This is also the subfoulder where the network input files are stored: 
                        #   "./SIMNAME/network_inputs_NETNAME.xlsx

NETNAME <- c("monarch") # Give a distinct name for each class as used in input files
# Order is important here we would index [[1]] = class 1 and [[2]] = class 2 in alpha and beta.

SAVENAME <- "CodeMonarchs_NewEqn" # Specifies which folder contains species specific data

PERTNAME <- "PERTmonarchs" # Specifies the name of the .r perturbation file

PATH_PERTURB <- "Perturb_Document_Code_Data" # Specifies which folder contains perturbation data

PERTSAVE <- "PERTdataMonarch" # Specifies the name of the .RData file where the perturbation data will be saved

RNW <- "MonarchAnalysisReviewJune2018" # Specifies .rnw file to use when producing PDF

ERR <- .01 # Error tolerance for convergence. 
# To test convergence, we compare total population of all classes in the current season
# to the matching season from the previous year.

seasons <- 7 # Number of seasons or steps in one annual cycle. 
# This must match number of spreadsheets in input files
num_nodes <- 4 # Number of nodes in the network
# This must match the number of initial conditions given in input files

tmax <- 301 # Maximum number of steps to take - assume non convergence if t=tmax

OUTPUTS <- TRUE # TRUE = Process final outputs, FALSE = Do not process just run the simulation.

delta <- .3 # Used in KR Calculation 1 = full node removal 0 = do nothing

SN_length <- matrix(c(6, 1, 1, 1, 1, 1, 1),1,seasons) # length of the seasons used in perturbation and KR calculation.

RUN_PERT <- TRUE #TRUE = Code will run perturbations with the below rates
                  # FALSE = Code will just run baseline simulation
PERT <- c(.9, .8, .7, .6, .5, 1) # survival rate perturbations desired if doing perturbations

## For debugging your model equations ##
SILENT <- TRUE # TRUE = Do not print data to console - silence outputs.
                # FALSE = Print population data and network function data to the Console for debugging.
SAVE_VAR <- TRUE

### Users should not need to interact with the code below ###

################
## SIMULATION ##
################

# Set the working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Set location of source code
netcode <- c("../NetworkCode1.2_NewEqn/")
# Clear the workspace reserving needed network input variables
base_variables <- c("seasons", "num_nodes", "NETNAME", "tmax", "SIMNAME", "ERR", "OUTPUTS", "SILENT","netcode","base_variables", "SAVE_VAR", "RUN_PERT", "delta", "SN_length", "PERT","SAVENAME","RNW","PERTNAME","PATH_PERTURB","PERTSAVE")
if(!exists("pert_variables")){pert_variables <- c("pert_variables")}
pert_variables <- c(pert_variables, ls(),"count", "p") 
base_variables <- c(base_variables, "pert_variables")
rm(list=setdiff(ls(), base_variables))
if(SAVE_VAR == TRUE){
  save_variables <- c("save_variables")
  base_variables <- c(base_variables, "save_variables")   }


### SET UP THE NETWORK(S) ###
source(paste(netcode,"NetworkSetup.R",sep=""))

### RUN THE BASELINE SIMULATION ###
print(paste("Running", SIMNAME, sep=" "))
source(paste(netcode,"NetworkSimulation.R",sep=""))

########################  
###  PROCESS OUTPUTS ###
########################

if (OUTPUTS == T){
  source(paste(netcode,"NetworkOutputs.R",sep=""))
}

########################
##  RUN PERTURBATIONS ##
########################
if (RUN_PERT == T){
  source(paste("../", SAVENAME, "/",PERTNAME,".R",sep = ""))
}


######################################
### SAVE THE DATA TO RUN .RNW FILE####
######################################

save.image(file=paste("../", SAVENAME,"/",PATH_PERTURB,"/DataFiles/", PERTSAVE, ".RData",sep="" ))


##########################
### COMPILE PDF FILE   ###
##########################

Sweave(paste("../",SAVENAME,"/",PATH_PERTURB,"/",RNW,".rnw",sep=""), encoding = "UTF-8")

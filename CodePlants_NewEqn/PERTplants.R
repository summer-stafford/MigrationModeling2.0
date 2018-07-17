## RUN the SURVIVAL PERTURBATIONS
## Run ElkSimulation.R before running Perturbations
# Set the working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# SAVE ORIGINAL VALUES:
# SAVE ORIGINAL VALUES:
BASECR <- CR
BASECRt <- CRt
BASECRs <- CRs
BASEWR <- WR
BASEWRt <- WRt
BASEWRs <- WRs
BASELAMBDAt <- LAMBDAt
BASEGAMMAt <- GAMMAt
BASEKR <- KR
BASEKRs <- KRs
BASEDS <- DS
BASEPop_TOTAL <- Pop_TOTAL
BASEN <- N
BASEtimestep <- timestep
BASEtotal_pop <- total_pop



print("Perturbing Survival Rates at each Node")
print(PERT)

PERT_CR <- list()
PERT_CRt <- list()
PERT_CRs <- list()
PERT_WR <- list()
PERT_WRt <- list()
PERT_WRs <- list()
PERT_LAMBDAt <- list()
PERT_GAMMAt <- list()
PERT_KR <- list()
PERT_KRs <- list()
PERT_DS <- list()
PERT_Pop_TOTAL <- list()
PERT_timestep <- list()
PERT_N <- list()
PERT_storerun <- list()
BASE_alpha <- alpha

if(SAVE_VAR==TRUE){
  BASE_save_f_update <- save_f_update
  BASE_save_p_update <- save_p_update
  BASE_save_s_update <- save_s_update
  BASE_save_CR <- save_CR 
  
  PERT_save_f_update <- list()
  PERT_save_p_update <- list()
  PERT_save_s_update <- list()
  PERT_save_CR <- list()
  save_variables <- c(save_variables, "BASE_save_f_update", "BASE_save_p_update", "BASE_save_s_update", "BASE_save_CR", "PERT_save_f_update", "PERT_save_p_update", "PERT_save_s_update", "PERT_save_CR")
}

pert_variables <- c(pert_variables, ls(),"count", "p")
count <- 0

for (p in 1:length(PERT)){
  for (num in 1:num_nodes){
    count <- count + 1
  # Update network parameters:
  # alpha[[class]][[season]]$S - contains baseline survival rates at all nodes
  # alpha[[class]][[season]]$smin - contains density dependent paramters for post harvest survival calculation
  # alpha[[class]][[season]]$smax - contains density dependent paramters for post harvest survival calculation
  # NOTE: post harvest survival is updated as: f_new[type,type] <- (smin+((smax)-(smin))/(1+exp(-Z))) so perturbing both smax and smin equally will perturb this survival rate.
  alpha <- BASE_alpha # reset alpha to baseline
  for (cl in 1:NUMNET){
    # Reset N-initial to the last (season 1) timestep in the simulation
    N_initial[[cl]][num] <- N[((num-1)*NUMNET+cl) ,timestep]
    # If commented out this resets N-initial to the original initial populations in the baseline.
    
    for (SN in 1:seasons){  # perturb all seasons of breeding and winter !!!!!!!!!!! THIS MUST BE CHANGED FOR EACH SPECIES !!!!!!!!!!!!!!
        alpha[[cl]][[SN]]$s_node[[num]] <- BASE_alpha[[cl]][[SN]]$s_node[[num]]*PERT[p]
        # CANT MAKE SURVIVAL RATES GREATER THAN 1
        if(alpha[[cl]][[SN]]$s_node[[num]]>1){ 
          print(paste("Survival rate for PERT: ", PERT[p], "on Node", num, "season", SN, "is",alpha[[cl]][[SN]]$s_node[[num]], ". Resetting to 1", sep=" " ))
          alpha[[cl]][[SN]]$s_node[[num]] <- 1 
          }
    }
  }
  
  # Make sure we are not starting the population at zero
  if(sum(N_initial[[1]][1:3])<1000){ N_initial[[1]][1:3] <- c(1000,1000,1000)}
  if(sum(N_initial[[2]][1:3])<1000){ N_initial[[2]][1:3] <- c(1000,1000,1000)}
  
  # RUN SIM
  print(paste("Running PERT", PERT[p], "on Node", num, sep=" "))
  PERT_storerun[[count]] <- paste("Running PERT", PERT[p], "on Node", num, sep=" ")
  source(paste(netcode,"NetworkSimulation.R",sep=""))
  source(paste(netcode,"NetworkOutputs.R",sep=""))
  
  # Save Results
  PERT_CR[[count]] <- CR
  PERT_CRt[[count]] <- CRt
  PERT_CRs[[count]] <- CRs
  PERT_WR[[count]] <- WR
  PERT_WRt[[count]] <- WRt
  PERT_WRs[[count]] <- WRs
  PERT_LAMBDAt[[count]] <- LAMBDAt
  PERT_GAMMAt[[count]] <- GAMMAt
  PERT_KR[[count]] <- KR
  PERT_KRs[[count]] <- KRs
  PERT_Pop_TOTAL[[count]] <- Pop_TOTAL
  PERT_timestep[[count]] <- timestep
  PERT_N[[count]] <- N
  PERT_save_f_update[[count]] <- save_f_update
  PERT_save_p_update[[count]] <- save_p_update
  PERT_save_s_update[[count]] <- save_s_update
  PERT_save_CR[[count]] <- save_CR 
  }
}
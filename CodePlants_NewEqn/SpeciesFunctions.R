#################################
##    MODEL FUNCTIONS          ##
#################################

### PLANT ###

############## CALCULATE f -> NODE DYNAMICS UPDATE #################

f_function <- function() {
  ## This function must output the species demographic and class transition matrix
  # TYPICAL VARIABLES USED (N,alpha,i,ind,t)
  # N - contains the population arriving at the node at all previous time steps - N["time_step"] ordered by class
  # NOTE if you specify "class_name" then it will look to that specific class's population each time
  # NOTE if you set "class_name"=[[i]] it will get the numbers for the current calculation's class
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # i - gives the class number as an integer - the class is are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then i 1 = seeds and i 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  
  # Initiaize basic update information
  f_new <- matrix(0,nrow=NUMNET, ncol=NUMNET)
  
  # Density dependent equations:
  
  ### RETRIEVE F UPDATE EQUATIONS ###
  Survi_seed <- alpha[[1]][[ind]]$s_node[node]
  Survi_plant <- alpha[[2]][[ind]]$s_node[node]
  trans <- alpha[[2]][[ind]]$t_node[node]
  repro <- alpha[[1]][[ind]]$r_node[node]
  p <- alpha[[2]][[ind]]$psi_node[node] #psi value at node
  POP_seeds <- N[node*2-1,t]
  POP_plants <- N[node*2,t]
  K <- alpha[[2]][[ind]]$K_node[node]
  psi <- p * exp(-(POP_plants+trans*Survi_seed*POP_seeds)/K) #density dependent transition equation
  
  for (type in 1:NUMNET) {
    ### UPDATE SEED POPULATION ###
    if (type==1){
      f_new[type,type] <- Survi_seed*(1 - trans)
      f_new[2,type]<- repro
    }
    ### UPDATE PLANT POPULATION ###
    if (type == 2){
      f_new[type,type] <- Survi_plant
      f_new[1,type] <- psi*trans*Survi_seed
    }
  }
  
  # Units of f_new should be total number of migrants leaving the node
  
  return(f_new)
}
########################################################################## 

############## CALCULATE p -> PATH TRANSITION RATES #################
p_function <- function(){
  # TYPICAL VARIABLES USED (N,f_update,alpha,beta,i,ind,t)
  # N - contains the population arriving at the nodes at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE this population count is before accounting for node dynamics
  # f_update - contains the population after node dyanamics have been applied for the current time step - f_update[["class_number"]]
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # beta - contains all information about path dynamics - beta[["class_number"]][["season]]$"variable_name"
  # NOTE edge transition rates are found in beta[["class_number"]][["season]]$p_ij
  # i - gives the class number as an integer - the class is are defined in NETNAME 
  # for example: if NETNAME <- c("seeds","plants") then i 1 = seeds and i 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  # NOTE if it is required to change any variable globally then use "old_value" <<- "new_value"
  
  p <- beta[[type]][[ind]]$p_ij 
  p_new <- p
  # p_new should be unitless and represent the edge transition probabilities for the given step
  
  return(p_new)
}
##########################################################################

############## CALCULATE s -> PATH SURVIVAL RATES #################
s_function <- function(){
  # TYPICAL VARIABLES USED (N,f_update,alpha,beta,i,ind,t)
  # N - contains the population arriving at the nodes at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE this population count is before accounting for node dynamics
  # f_update - contains the population after node dyanamics have been applied for the current time step - f_update[["class_number"]]
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # beta - contains all information about path dynamics - beta[["class_number"]][["season]]$"variable_name
  # NOTE edge survival rates are found in beta[["class_number"]][["season]]$s_ij
  # i - gives the class number as an integer - the class is are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then i 1 = seeds and i 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  # NOTE if it is required to change any variable globally then use "old_value" <<- "new_value"
  
  # Edge survival where k=hunting mortality
  s <- beta[[type]][[ind]]$s_ij
  
  s_new <- s
  # s_new should be unitless and represnt the edge survival for the given step
  return(s_new)
}
##########################################################################

#################################
##    MODEL FUNCTIONS          ##
#################################

#################################
##             ELK             ##
#################################

############## CALCULATE f -> NODE DYNAMICS UPDATE #################

f_function <- function() {
  # TYPICAL VARIABLES USED (N,alpha,type,ind,t)
  # N - contains the population arriving at the node at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE if you specify "class_name" the it will look to that specific class's population each time
  # NOTE if you set "class_name"=[[type]] it will get the numbers for the current calculation's class
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  
  
  # Create f_new matrix to store new population values
  f_new <- matrix(0, nrow= NUMNET, ncol= NUMNET)
  
  # Density dependent equations:
  SI_A <- alpha[[2]][[ind]]$s0[[node]] # Baseline adult survival rate type=2
  SI_C <- alpha[[1]][[ind]]$s0[[node]] # Baseline calf survival rate type=1
  K <- alpha[[1]][[ind]]$K[[node]]+alpha[[2]][[ind]]$K[[node]] # Total carying capacity
  Na <- N[node * 2, t]# Number of adults pre
  Nc <- N[node * 2 -1,t] # Number of calves pre
  R <- alpha[[2]][[ind]]$r[node] # Adult reproductive rate
  Survi_juvenile <- SI_C*exp(1-(Na+Nc)/K) #Density dependent survival rate in the summer and fall for nodes 1&3
  Survi_adult <- sqrt(SI_A*exp(-0.219*((Na+Nc)/K)^3.77))
  
  #Breeding nodes
  if (ind == 2){ #Breeding season ind = 2 because of the way alpha is set up
    for (type in 1:NUMNET){ 
      ################# UPDATE ADULT POPULATION ####################
      if (type == 2){
        f_new[type,type] <- Survi_adult
        f_new[1,type] <- Survi_juvenile #Transition rate is 1 in the breeding season, so all surviving juveniles transition
      }
      ################ UPDATE JUVENILE POPULATION #################
      if (type == 1){
        f_new[type,type] <- 0 #there are no surviving calves that do not trasition
        f_new[2,type] <- 0.6*R*Survi_adult
      }
    }
  }
  
  #Nonbreeding nodes
  if (ind == 1){
    for (type in 1:NUMNET){
     if (type == 2){
       f_new[type,type] <- Survi_adult
       f_new[1,type] <- 0 #Transition rate in season two is zero
     }
      if (type == 1){
        f_new[type, type] <- SI_C  #constant during the winter
        f_new[2,type] <- 0 #non breeding season
      }
    }
  }
  
  return(f_new)
}
########################################################################## 

############## CALCULATE p -> PATH TRANSITION RATES #################
p_function <- function(){
  # TYPICAL VARIABLES USED (N,f_update,alpha,beta,type,ind,t)
  # N - contains the population arriving at the nodes at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE this population count is before accounting for node dynamics
  # f_update - contains the population after node dyanamics have been applied for the current time step - f_update[["class_number"]]
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # beta - contains all information about path dynamics - beta[["class_number"]][["season]]$"variable_name"
  # NOTE edge transition rates are found in beta[["class_number"]][["season]]$p_ij
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # for example: if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  # NOTE if it is required to change any variable globally then use "old_value" <<- "new_value"
  
  p <- beta[[type]][[ind]]$p_ij
  n <- length(p)
  resident <- matrix(0,n,n)
  
  if(ind == 1){ # winter
    Npre3 <- N[type+NUMNET*2,t] #type +NUMNET*(node-1) where node = 3
    Npre2 <- N[type+NUMNET,t]
    if(t==1){
    
      resident[3,3] <- Npre3-0.13/0.87*Npre2
    }
    else{
      resident[3,3] <- M[[t-1]][[type]][3,3]
    }
    
    if(Npre3 > 0){
      p[3,3] <- resident[3,3]/Npre3
      p[3,1] <- 1 - resident[3,3]/Npre3
    }
    else{
      p[3,3] <- 0
      p[3,1] <- 0
    }
  }
  
  p_new <- p
  
  return(p_new)
  # p_new should be unitless and represent the edge transition probabilities for the given step
  
}
##########################################################################

############## CALCULATE s -> PATH SURVIVAL RATES #################
s_function <- function(){
  # TYPICAL VARIABLES USED (N,f_update,alpha,beta,type,ind,t)
  # N - contains the population arriving at the nodes at all previous time steps - N[["time_step"]]$"class_name"
  # NOTE this population count is before accounting for node dynamics
  # f_update - contains the population after node dyanamics have been applied for the current time step - f_update[["class_number"]]
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # beta - contains all information about path dynamics - beta[["class_number"]][["season]]$"variable_name
  # NOTE edge survival rates are found in beta[["class_number"]][["season]]$s_ij
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  # NOTE if it is required to change any variable globally then use "old_value" <<- "new_value"
  
  s <- beta[[type]][[ind]]$s_ij
  
  s_new <- s
  # s_new should be unitless and represent the edge survival for the given step
  
  return(s_new)
}
##########################################################################

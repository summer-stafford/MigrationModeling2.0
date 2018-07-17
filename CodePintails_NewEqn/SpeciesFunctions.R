#################################
##    MODEL FUNCTIONS          ##
#################################

### PINTAILS ###

############## CALCULATE f -> NODE DYNAMICS UPDATE #################

f_function <- function() {
  ## This function must output the species demographic and class transition matrix
  # TYPICAL VARIABLES USED (N,alpha,type,ind,t)
  # N - contains the population arriving at the node at all previous time steps - N["time_step"] ordered by class
  # NOTE if you specify "class_name" the it will look to that specific class's population each time
  # NOTE if you set "class_name"=[[type]] it will get the numbers for the current calculation's class
  # alpha - contains all information about node dynamics - alpha[["class_number"]][["season"]]$"variable_name"
  # type - gives the class number as an integer - the class types are defined in NETNAME 
  # NOTE if NETNAME <- c("seeds","plants") then type 1 = seeds and type 2 = plants
  # ind - gives the season number as an integer
  # t - gives the time step
  
  # Initiaize basic update information
  f_new <- matrix(0,nrow=NUMNET, ncol=NUMNET)
  
  ##########ERASE THIS###############
  # # Initiaize basic update information
  # for (type in 1:NUMNET){
  #   f_new[type,type] <- alpha[[type]][[ind]]$S[[node]]
  # }
  
  
  
  # Breeding Nodes
  if(ind==1){
    for (type in 1:NUMNET){
      
      ########ADD THIS##############
      if(type==1 || type==2){ # If we are updating the Adult population
        f_new[type,type] <- alpha[[type]][[ind]]$S[[node]]
      }
      if(type==3 || type==4){ # If we are updating the Juvenile population new births
        
        ########CHANGE POP DEFINITION##############
        #POP <- sum(N[seq(1,length(N), by=num_nodes)],t)+sum(N[seq(2,length(N), by=num_nodes)],t) # Number of existing Males and Females for density dependence
        POP <- sum(N[seq(NUMNET*(node-1)+1,node*NUMNET),t])
        
        Survi_F <- alpha[[1]][[ind]]$S[[node]] # Female survival Rate
        
        ponds <- alpha[[type]][[ind]]$P[[node]]
        alpha0 <- alpha[[type]][[ind]]$a_0[[node]]
        alpha1 <- alpha[[type]][[ind]]$a_1[[node]]
        alpha2 <- alpha[[type]][[ind]]$a_2[[node]]
        
        R <- exp(alpha0+alpha1*POP+alpha2*ponds) # Density Dependent Reproductive Rate
        
        f_new[1,type] <- R*Survi_F
        
      } 
    }
  }
  
  
  # Winter Nodes
  # Density dependent node survival rates only apply to the post harvest winter-spring season
  if(ind==2){
    for (type in 1:NUMNET){
      beta0 <- alpha[[type]][[ind]]$b_0[[node]]
      beta1 <- alpha[[type]][[ind]]$b_1[[node]]
      smax <- alpha[[type]][[ind]]$S_max[[node]]
      smin <- alpha[[type]][[ind]]$S_min[[node]]
      ########CHANGE POP DEFINITION##############
      #POP <- sum(N)
      POP <- sum(N[seq(NUMNET*(node-1)+1,node*NUMNET),t])
      
      if(type==1 || type==2){
        Z <- beta0+beta1*POP 
        f_new[type,type] <- (smin+((smax)-(smin))/(1+exp(-Z))) # Density Dependent Survival
      }
      if(type==3){
        Survi_J <- (smin+((smax)-(smin))/(1+exp(-Z)))
        f_new[type,1] <- Survi_J # Transition the Juveniles to Adult Females
      }
      if(type==4){
        Survi_J <- (smin+((smax)-(smin))/(1+exp(-Z)))
        f_new[type,2] <- Survi_J # Transition the Juveniles to Adult Males
      }
    }
  }
  ########ADD THIS##############
  #Stopover nodes
  if(ind==3){
    for (type in 1:NUMNET){
      if(type==1 || type==2){ # If we are updating the Adult population
        f_new[type,type] <- alpha[[type]][[ind]]$S[[node]]
      }
    }
  }
  # f_new should be a NUMNET X NUMNET matrix containing demographic rates
  
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
  
  p <- beta[[type]][[ind]]$p_ij  #p_edge[[type]][[ind]]
  
  # density dependent transitions only apply to the spring stopover season
  if(ind==3){
    delta0 <- alpha[[type]][[ind]]$Delta_0[2]
    delta1 <- alpha[[type]][[ind]]$Delta_1[2]
    delta2 <- alpha[[type]][[ind]]$Delta_2[2]
    
    # Number of migrants who survived PR
    PRPOP <- sum(f_update[[2]] %*%  N[seq(NUMNET+1,2*NUMNET,by=1),t])
    ponds <- alpha[[type]][[ind]]$P[2]
    
    Y <- delta0 + delta1*PRPOP+delta2*ponds
    
    # Original transition probabilities from PR
    psim <- alpha[[type]][[ind]]$psi_max[2]
    newpsi_leave <- (psim)*1/(1+exp(-Y))
    
    # Probabilities of leaving PR for AK or NU
    psi21 <- beta[[type]][[ind]]$psi_ij[2,1]  
    psi23 <- beta[[type]][[ind]]$psi_ij[2,3]
    
    # Update the transition probabilities
    p[2,1] <- newpsi_leave*psi21
    p[2,3] <- newpsi_leave*psi23
    p[2,2] <- 1-newpsi_leave
    
  }
  
  p_new <- p
  # p_new should be unitless and represent the edge transition probabilities for the given step
  
  return(p_new)
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
  
  # Edge survival where k=hunting mortality
  s <- beta[[type]][[ind]]$s_ij*(1-beta[[type]][[ind]]$kappa_ij)  #s_edge[[type]][[ind]]
  
  s_new <- s
  # s_new should be unitless and represnt the edge survival for the given step
  return(s_new)
}
##########################################################################
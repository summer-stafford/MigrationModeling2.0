#######################################
##           SIMULATIONS             ##
#######################################

## This code runs the network simulation based on the general network model
## It first loads the user defined Newtork Functions: f_(i,t), p_(ij,t), and s_(ij,t)

SAVE_VAR <- TRUE
## USERS SHOULD NOT NEED TO INTERACT WITH THIS CODE
AMATRIX_CR <- list()

source("SpeciesFunctions.R")

f_update <- list()
p_update <- list()
s_update <- list()

if(SAVE_VAR == TRUE){
  save_f_update <- list()
  save_p_update <- list()
  save_s_update <- list()
  save_CR <- list()
  save_variables <- c(save_variables, "save_f_update", "save_p_update", "save_s_update", "save_CR")
  SEASONNAMES <- c()
  for (i in 1:seasons){
    SEASONNAMES <- c(SEASONNAMES, paste("season", i))
  }
  crtime <- 0
}

## CREATE M MATRIX ##
M <- list()
## INITIALIZE NODE POPULATION
N <- matrix(0,nrow=num_nodes*NUMNET, ncol=1) #Set initial population after movement to each node (before reproduction/survival at node)

for (j in 1:NUMNET){
  for (i in 1:num_nodes){
  N[j + NUMNET*(i-1)] = N_initial[[j]][i]
  }
}

### INITIALIZE OTHER POPULATION DATA
total_pop <- matrix(0,nrow=NUMNET, ncol=1) # total network population - sum of class population at each node
rownames(total_pop) <- NETNAME
proportional_pop <- matrix(0,nrow=num_nodes, ncol=1) # proportion of 
## rownames(proportional_pop) <- named after nodes

for (i in 1:NUMNET){
  total_pop[i] <- sum(N[seq(i,length(N), by=NUMNET)]) 
}
rm(i)

if(SILENT==F){
  print("Initial Population")
  print(total_pop)
}

### INITIALIZE ERROR ###
errorstop <- 0
ERRPOP_OLD <- matrix(100,1,seasons)  

### INITIALIZE TIME STEPS ###
t <- 0


### START SIMULATION ###
### LOOP THROUGH TIME ###
while (errorstop == 0){
  t <- t+1
  if(SILENT==F){
    print("Time Step")
    print(t)
  }
  
  ind <- ((t-1)%%seasons) +1 # season number
  if(SILENT==F){
    print("Season Number")
    print(ind)
  }
  
  ## MODEL FUNCTIONS ##
  for (i in 1:num_nodes){
    node <- i
    f_update[[i]] <- f_function()   # f_function(N,alpha,i,ind,t) # Node update demographics dimensions = class X class
  }
  #print(f_update)
  for (i in 1:NUMNET){
    type <- i
    p_update[[i]] <- p_function()    # p_function(N,f_update,alpha,p_edge,beta,i,ind,t) # Path transition probability dimensions = num_node X num_node
    s_update[[i]] <- s_function()   # s_function(N,f_update,alpha,s_edge,beta,i,ind,t) # Path survival probability dimensions = num_node X num_node
  }
  rm(i)
  if(SAVE_VAR == TRUE){
    save_f_update[[t]] <- f_update
    save_p_update[[t]] <- p_update
    save_s_update[[t]] <- s_update
  }
  
  if(SILENT==F){
    print("f_(i,t)=")
    print(f_update)
    print("p_(ij,t)=")
    print(p_update)
    print("s_(ij,t)=")
    print(s_update)
  }
  
  ## MODEL EQUATION ##
  ## Create matrix columns for next time step
  N <- cbind(N,matrix(0,nrow=NUMNET*num_nodes,ncol=1))
  
  ## MATRICEES NEEDED FOR UPDATE ##
  ## F Block Diagonal ##
  FBLOCK <- matrix(0,nrow = NUMNET*num_nodes, ncol = NUMNET*num_nodes)
  for (i in 1:num_nodes){
    E <- matrix(0,nrow = num_nodes, ncol = num_nodes)
    E[i,i] <- 1
    FBLOCK <- FBLOCK + kronecker(E, t(matrix(f_update[[i]],nrow=NUMNET,ncol=NUMNET)))
  }
  rm(i,E)
  ## Q Block Diagonal ##
  
  QBLOCK <- matrix(0,nrow = NUMNET*num_nodes, ncol = NUMNET*num_nodes)
  for (i in 1:NUMNET){
    E <- matrix(0,nrow = NUMNET,ncol = NUMNET)
    E[i,i] <- 1
    QBLOCK <- QBLOCK + kronecker(t(matrix(unlist(s_update[[i]]),nrow = num_nodes,ncol = num_nodes) * matrix(unlist(p_update[[i]]),nrow = num_nodes,ncol = num_nodes)),E)
  }
  rm(i,E)

  AMATRIX <- QBLOCK %*% FBLOCK
  AMATRIX_CR[[t]] <- AMATRIX
  
    if(SILENT==F){
    print("FBLOCK")
    print(FBLOCK)
    print("QBLOCK")
    print(QBLOCK)
    print("A")
    print(AMATRIX)
  }
  
  ######################## UPDATE MIGRANTS LIST#######################################
  M[[t]] <- list()
  for(j in 1:NUMNET){
    ### CREATE MATRIX AT EACH TIME STEP FOR BOTH CALVES AND ADULTS
    M[[t]][[j]] <- matrix(0, nrow= num_nodes, ncol = num_nodes)
    
    matrix_c <- matrix(0, nrow= num_nodes * NUMNET, ncol=num_nodes)
    for (i in 1:num_nodes){
      E <- matrix(0,nrow = num_nodes, ncol = num_nodes)
      E[i,i] <- 1
      e <- matrix(0,nrow = num_nodes, ncol = 1)
      e[i] <- 1
      matrix_c <- matrix_c + (kronecker(E, t(matrix(f_update[[i]],nrow=NUMNET,ncol=NUMNET))))%*%N[,t]%*%t(e)
    }
    rm(E,i,e)
  
    matrix_b <- matrix(0, nrow= num_nodes * NUMNET, ncol= num_nodes * NUMNET)
    E <- matrix(0,nrow = NUMNET,ncol = NUMNET)
    E[j,j] <- 1
    matrix_b <- kronecker(t(matrix(unlist(s_update[[j]]),nrow = num_nodes,ncol = num_nodes) * matrix(unlist(p_update[[j]]),nrow = num_nodes,ncol = num_nodes)),E)
    rm(E)
  
    matrix_a <- matrix(0, nrow = num_nodes, ncol = num_nodes * NUMNET)
    matrix_ident <- diag(num_nodes)
    vector_class <- matrix(1,nrow = NUMNET, ncol=1)
    matrix_a <- kronecker(matrix_ident, t(vector_class))
    rm(matrix_ident,vector_class)
  
    M[[t]][[j]] <- t(matrix_a %*% matrix_b %*% matrix_c)
    #print(matrix_a)
    #print(matrix_b)
    #print(matrix_c)
    rm(matrix_a,matrix_b,matrix_c)
  }
  ## Calculate next time step
  N[,t+1] <- AMATRIX %*% N[,t]
  
  
  ## UPDATE TOTAL POPULATION ##
  total_pop <- cbind(total_pop,matrix(0,nrow=NUMNET,ncol=1))
  ## Set total population data
  for (i in 1:NUMNET){
    curpop <- N[,t+1]
    total_pop[i,t+1] <- sum(curpop[seq(i,length(curpop), by=NUMNET)])
  }
  
  ## Cr as a function of time ##
  if(SAVE_VAR == TRUE){
    if(t >= 2*seasons){ # Need two seasons of data to calculate full year CR values
      if(ind==1){ # Reach breeding season, go back and calculate last year's CR values
        crtime <- crtime + 1
        startCRtime <- t - 1 - seasons
        source(paste(netcode,"CRKR_ONE.R",sep=""))
        save_CR[[crtime]] <- CR
        colnames(save_CR[[crtime]]) <- SEASONNAMES
        NAMESofROW <- c()
       for (i in 1:num_nodes){
          for (k in 1:NUMNET){
           NAMESofROW <- c(NAMESofROW, paste(colnames(beta$adult_female[[1]]$p_ij)[i], NETNAME[k]))
          }
        }
        NAMESofROW <- gsub("_","",NAMESofROW)
        row.names(save_CR[[crtime]]) <- NAMESofROW
      }
    }
  }
  
  if(SILENT==F){
    print(paste("Total Population in season", ind, "for time step", t))
    print(total_pop[,t+1])
  }
  
  if(SILENT==F){
    print(paste("Node Population in season", ind, "for time step", t))
    print(N[,t+1])
  }
  
  ### TEST FOR BLOW-UP OR INFINITE NUMBERS ###
  if(any(is.finite(total_pop[,t+1]))=="FALSE"){
    stop("\n NaN found in population data: \n *** This means an infinite population has been reached in at least one node. \n *** Check divide by zero or missing data in model parameters. \n *** Check density dependent equations. \n *** This is likely caused by SpeciesFunctions.R")}
  
  
  ### CALCULATE THE ERROR AT EACH SEASON ###
  if(t >= seasons){
    sum_new <- sum(total_pop[,t+1])
    sum_old <- sum(total_pop[,t+1-seasons])
    ERRPOP_OLD[1,ind] <- abs(sum_new - sum_old)
    ## Only stop if the total error is less than the allowable error across all seasons
    if(all(ERRPOP_OLD < ERR)){errorstop <- 1}
    if(t+1 >= tmax){
      errorstop <- 1 
      cat("\n ************ \n\n The simulation did not converge within the maximum time allowed \n\n ************ \n")}
  }
}

# We give results for the total population at the start of season 1
# Check at which season the simulation stopped
timestep <- t+1 - ind%%seasons

# Clear the workspace reserving needed network input variables and base variables and simulation variables
simulation_variables <- c("timestep","total_pop","N","t","ind","f_update","p_update","s_update","simulation_variables", "AMATRIX", "AMATRIX_CR")
rm(list=setdiff(ls(),c(network_variables,base_variables,simulation_variables,save_variables,pert_variables)))
#need to add in pert_variables and save_variables


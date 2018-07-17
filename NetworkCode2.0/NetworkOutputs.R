#################################
##        BASIC OUTPUTS        ##
#################################

## This code produces basic outputs for the Network simulation.
## 1. A graph of Population over time for each class - to show each population reaching a steady state.
## 2. A table of values showing the steady state total network populations for each season at steady state.
## 3. A single season graph of the steady state total network populations for each season at steady state.

## USERS SHOULD NOT NEED TO INTERACT WITH THIS CODE
options(digits=4)
options(scipen=999)

# CALCULATE CR
startCRtime = timestep - 1 - seasons 
source(paste(netcode,"CRKR_ONE.R",sep=""))

# # CALCULATE Pop_TOTAL  --> maybe move this to CRKR_ONE.R
 Pop_TOTAL <- matrix(0,nrow = num_nodes, ncol = seasons)
 for (j in 1:seasons){
   for (i in 1:num_nodes){
     TIME <- timestep - 1 - 2*seasons+j
     Pop_TOTAL[i,j] <-sum(N[seq((i-1)*NUMNET+1,i*NUMNET),TIME]) 
   }
 }

#### Plot total populations over time  ####

# Only plot on the first run, not during perturbation experiments
if(!exists("PERT")){
  
platform <- .Platform$OS.type
if(platform == "unix") {X11()
} else{
  if(platform == "windows"){windows()
  } else{quartz()}}

time_steps <- 0:timestep/seasons
STEPwidth <- seasons
STEPstart <- 1
graphbot <- 0 
graphtop <- max(total_pop[seq(1,timestep,by=seasons)])
dottype <- data.frame(matrix(0,1,NUMNET))
dottype[1] <- 20

plot(c(0),c(0),type="l", 
     main=paste("Total Network Population vs. Time \n Single Annual Survey (beginning of season 1)"),
     ylim = c(graphbot,graphtop), xlim = c(min(time_steps),max(time_steps)),
     xlab="Years",ylab="Total Population")

xplot <- time_steps[seq(STEPstart, length(time_steps),by=STEPwidth)]

for(i in 1:NUMNET){
  par(new=TRUE)
  yplot <- total_pop[i,seq(STEPstart, length(time_steps),by=STEPwidth)]
  dot <- as.double(dottype[i])
  plot(xplot,yplot,type="o", lty=1, pch=dot,ylim = c(graphbot,graphtop), 
       xlim = c(min(time_steps),max(time_steps)), axes="FALSE", xlab = "", ylab = "")
  dottype[i+1] <- dottype[i] + 1
}
rm(i)

legend('right', NETNAME , lty=1, pch=dottype, bty='n', cex=.75)

}

# #### Store data as .csv for steady state annual cycle population numbers  ####  
# pop_output <- data.frame(matrix(0,seasons,NUMNET))
# for(i in 1:NUMNET){
#   pop_output[,i] <- data.frame(total_pop[seq(timestep-seasons+1,timestep, by=1)],row.names=NULL) 
# }
# rm(i)
# 
# colnames(pop_output) <- NETNAME
# for(i in 1:seasons){
#   rownames(pop_output)[i] <- paste("season",i)
#   }
# 
# 
# print("Total Equilibrium Population - for each season")
# print(pop_output)
# rm(i)
# 
# 
# pop_dist <- matrix(0,num_nodes,seasons)
# for(i in 1:seasons){
#   node_pop <- rowSums(matrix(unlist(N[[timestep-(seasons-i+1)]]), nrow=num_nodes))
#   pop_dist[,i] <- node_pop/sum(node_pop)
# }
# colnames(pop_dist) <- rownames(pop_output)
# rownames(pop_dist) <- colnames(M[[1]][[1]])
# cat("\n\n Node Equilibrium Population Distribution - for each season\n")
# print(pop_dist)
# 
# 
# 
# path_dist <- list()
# for(i in 1:seasons){
#   path_pop <- matrix(0,num_nodes,num_nodes)
#   for(j in 1:NUMNET){
#   path_pop <- path_pop + M[[timestep-(seasons-i+1)]][[j]]
#   }
#   path_dist[[i]] <- path_pop/sum(rowSums(path_pop))
#   rownames(path_dist[[i]]) <- colnames(path_dist[[i]])
# }
# cat("\n\n Path Equilibrium Population Distrubition - for each season\n")
# print(path_dist)
# 
# 
# #### Plot steady state annual cycle population numbers ####
# platform <- .Platform$OS.type
# if(platform == "unix"){X11()
# } else{
#   if(platform == "windows"){windows()}
#   else{quartz()}}
# 
# time_steps <- 1:seasons
# graphbot <- 0 
# graphtop <- max(pop_output)
# rm(dottype)
# dottype <- data.frame(matrix(0,1,NUMNET))
# dottype[1] <- 20
# 
# plot(c(0),c(0),type="l", ylim = c(graphbot,graphtop), xlim = c(min(time_steps),max(time_steps)),
#      xaxt = "n", main=paste("Class Population at Steady State \n Over One Annual Cycle"),
#      xlab="Season Number",ylab="Population")
# 
# axis(side = 1, at = time_steps)
# for(i in 1:NUMNET){
#   par(new=TRUE)
#   dot <- as.double(dottype[i])
#   plot(time_steps,pop_output[,i],type="o", lty=1, pch=dot, ylim = c(graphbot,graphtop), 
#        xlim = c(min(time_steps),max(time_steps)), axes = FALSE, xlab = "", ylab = "")
#   dottype[i+1] <- dottype[i] + 1
# }
# rm(i)
# legend('right', NETNAME , lty=1, pch=dottype, bty='n', cex=.75)
# 
# 
# #write.table(pop_output,file=paste(SIMNAME,"/SteadyStatePopulation.csv", sep=""), col.names=FALSE)
# 
# sink(file=paste(SIMNAME,"/SteadyStatePopulation.csv", sep=""))
# cat('Equilibrium Population')
# write.csv(pop_output)
# cat('___________________________')
# cat("\n")
# cat("\n")
# 
# cat('Equilibrium Population Distribution at Nodes')
# write.csv(pop_dist)
# cat('___________________________')
# cat("\n")
# cat("\n")
# 
# for(i in 1:seasons){
# cat('Equilibrium Population Distribution along Paths:   Season', i)
# write.csv(path_dist[[i]])
# cat('___________________________') 
# cat("\n")
# cat("\n")
# }
# 
# sink()

# Clear the workspace reserving needed network input variables and base variables and simulation variables
output_variables <- c("Pop_TOTAL", "SEASONNAMES", "CR" , "output_variables")
rm(list=setdiff(ls(),c(network_variables,base_variables,pert_variables,save_variables,simulation_variables,output_variables)))



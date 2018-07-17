
SEASONNAMES <- c()
for (i in 1:seasons){
  SEASONNAMES <- c(SEASONNAMES, paste("season", i))
}


## Calculate Cr values
POPnode <- matrix(0,nrow = num_nodes, ncol = seasons)
CR <- matrix(0,nrow = nrow(AMATRIX), ncol = seasons) #Cr values for each class/node/season
CRweight <- matrix(0,nrow = nrow(AMATRIX), ncol = seasons) #Cr*Nr values for each class/node/season 
ONES <- matrix(1,nrow = nrow(AMATRIX), ncol = 1)
CRt <- matrix(0,nrow = num_nodes, ncol = seasons) #Population weighted CR values for each node/season
#Choose each season as a focal season
for (SN in 1:seasons){
  CRtemp <- diag(nrow(AMATRIX))
  taustart <- startCRtime-seasons+SN  ## Calculation based on startCRtime
  tauend <- taustart+seasons-1
  for (i in taustart:tauend){
    CRtemp <- CRtemp %*% t(matrix(AMATRIX_CR[[i]], nrow = nrow(AMATRIX)))
  }
  CR[,SN] <- CRtemp %*% ONES
  CRweight[,SN] <- CR[,SN] * N[,taustart]
  
  for (j in 1:num_nodes){
    POPnode[j,SN] <- sum(N[((j-1)*NUMNET+1):(j*NUMNET),taustart])
    if(POPnode[j,SN]>0){
      CRt[j,SN] <- sum(CRweight[((j-1)*NUMNET+1):(j*NUMNET),SN])/POPnode[j,SN] # class based population weighted average within nodes
    }
  }
}

## Seasonal Population Weighted CR values
CRs <- matrix(0,nrow=num_nodes,ncol=1)
CRs <- rowSums(CRt*POPnode)/rowSums(POPnode)

if(SILENT == F){
  print("CR values vec{C}_t")
  print(CR)
  print(" ")
  print("Population weighted CR values (vec{C}_t * N)/Ntot")
  print(CRt)}
colnames(CR) <- SEASONNAMES
NAMESofROW <- c()
for (i in 1:num_nodes){
  for (j in 1:NUMNET){
    NAMESofROW <- c(NAMESofROW, paste(colnames(beta$adult_female[[1]]$p_ij)[i], NETNAME[j]))
  }
}
NAMESofROW <- gsub("_","",NAMESofROW)
row.names(CR) <- NAMESofROW
colnames(CRt) <- colnames(CR)
rownames(CRt) <- colnames(beta$adult_female[[1]]$p_ij)




## Population Proportion and Network Growth Rate
LAMBDAt <- matrix(0, 1, seasons)
WR <- matrix(0, NUMNET*num_nodes, seasons) #Contains population porportion for each class at each node.
WRt <- matrix(0,num_nodes,seasons) #Contains the population proportiona for each node, summing the classes
WRs <- matrix(0,num_nodes,1) #Contains the network wide average WR values - divide by number of seasons
for (SN in 1:seasons){
  TIME <- startCRtime-seasons+SN
  WR[,SN] <- N[,TIME]/sum(N[,TIME])
  TIME <- TIME + 1
}
colnames(WR) <- colnames(CR)
rownames(WR) <- rownames(CR)

for (i in 1:num_nodes){
  if(NUMNET==1){
    WRt[i,] <- WR[i,]
  }
  else{
    WRt[i,] <- colSums(WR[((i-1)*NUMNET+1):(i*NUMNET),])
  }
}

WRs <- rowSums(WRt)/3

for (SN in 1:seasons){
  LAMBDAt[SN] <- t(WR[,SN])%*%CR[,SN]
}
rownames(LAMBDAt)<-c("network growth rate")
colnames(LAMBDAt)<-colnames(CR)




## Kr and network growthrate in absense of r, GAMMA

KR <- matrix(0,num_nodes,seasons)
KRs <- matrix(0,num_nodes,1)
GAMMAt <- matrix(0,num_nodes,seasons)

Inc <- diag(1, nrow(AMATRIX), ncol(AMATRIX))
for (j in 1:num_nodes){
  Drr <- matrix(0, nrow(AMATRIX), ncol(AMATRIX))
  Enc <- matrix(0, nrow(AMATRIX), ncol(AMATRIX))
  xstart <- j*NUMNET-NUMNET+1
  xend <- j*NUMNET
  for (x in xstart:xend){
    TEMP <- matrix(0,nrow(AMATRIX), ncol(AMATRIX))
    TEMP[x,x] <- delta # Delta defined functional KR value NEED -- added season length weight
    Enc <- Enc + TEMP
  }
  Drr <- Inc - Enc

  
  
  for (SN in 1:seasons){
    CRnewtemp <- diag(nrow(AMATRIX))
    taustart <- startCRtime-seasons+SN  ## Calculation based on startCRtime
    tauend <- taustart+seasons-1
    for (i in taustart:tauend){
      CRnewtemp <- CRnewtemp %*% t(matrix(AMATRIX_CR[[i]], nrow = nrow(AMATRIX)))%*%Drr
    }
    GAMMAt[j,SN] <- WR[,SN]%*%Drr%*%CRnewtemp %*% ONES
  }
}

for (j in 1:num_nodes){
  KR[j,] <- LAMBDAt - GAMMAt[j,]
}
rownames(KR) <- colnames(beta$adult_female[[1]]$p_ij)
colnames(KR) <- colnames(CR)

POPsums <-colSums(POPnode)
TEMP <- matrix(0,num_nodes,seasons)
for (i in 1:seasons){
  TEMP[,i] <- matrix(POPsums[i],1,num_nodes)
}
KRs <- rowSums(KR*TEMP)/sum(POPsums)


## Ds - Proportional Dependence
DStemp <- matrix(0,num_nodes,seasons)
DS <- matrix(0,num_nodes,1)

DSdenom <- colSums(CRweight)
k <- seasons
# #CRweight contains Cr*Pt,i - weighted contributions.
for(i in 1:num_nodes){
  start <- (i-1)*NUMNET+1
  end <- i*NUMNET
  for(j in 1:seasons){
    DStemp[i,j] <- sum(CRweight[start:end,j])/DSdenom[j]  
  }
}
DS[,] <- 1/k * rowSums(DStemp)
rownames(DS) <- colnames(beta$adult_female[[1]]$p_ij)
colnames(DS) <- c("DS")


pert_variables <- c(pert_variables, "CR", "CRt", "CRs", "CRweight", "WR", "WRt", "WRs", "LAMBDAt","GAMMAt","KR", "KRs", "DS")
##############################
# simulation-calculation.R
##############################

library(divDyn)
library(paleotree)
library(foreach)
library(doParallel)

data(stages)
neg_stages <- -1 * stages$mid


#########################
## Run simulation
#########################

cores <- detectCores()
cl <- makeCluster(cores-1)
registerDoParallel(cl)

simulations <- foreach(i=1:500, .packages="paleotree") %dopar% {
  
simFossilRecord(p = 0.25, q = 0.25, r = 0.75, startTaxa = 200, 
                       totalTime = 66, nTotalTaxa = c(0,Inf), 
                       nExtant = c(0,Inf), nSamp = c(0,Inf))
  
}

stopCluster(cl)
save(simulations, file="simulations.RData")
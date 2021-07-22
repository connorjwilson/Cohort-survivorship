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

##########################
## Create FAD-LAD Matrix
##########################

sim_ranges <- foreach(i=1:500, .packages="paleotree") %dopar% {
  
  as.data.frame(fossilRecord2fossilRanges(simulations[[i]], ranges.only=TRUE))

}

sim_ranges_notna <- foreach(i=1:500) %dopar% {
 sim_ranges[[i]][!is.na(sim_ranges[[i]]$FAD),]
}


data(stages)
classify <- function (input) {
  
  output <- input
  
  for (i in 1:nrow(input)) {
    output[i,1] <- which((input[i,1]<=stages$bottom & input[i,1]>=stages$top)==TRUE)
    output[i,2] <- which((input[i,2]<=stages$bottom & input[i,2]>=stages$top)==TRUE)
  }
  
  return(output)
  
}

sim_ranges_stg <- foreach(i=1:500) %dopar% {
  classify(sim_ranges_notna[[i]])
}

for (i in 1:500) {
  sim_ranges_stg[[i]]$duration <- (sim_ranges_stg[[i]]$LAD - sim_ranges_stg[[i]]$FAD)

}

sim_surv  <- foreach(i=1:500, .packages="divDyn") %dopar% {
  survivors(fl=sim_ranges_stg[[i]])
}

sim_cohort <- list()
for (i in 1:500) {
  sim_cohort[[i]] <- sim_surv[[i]][82:95,82:95]
}



sim_fits <- foreach(i=1:500) %dopar% {
  fits <- list()
  for (j in 1:14) {
    fits[[j]] <- lm(log(sim_cohort[[i]][j:max(which(sim_cohort[[i]][,j]>0)),j])
                    ~ neg_stages[(81+j):(81+max(which(sim_cohort[[i]][,j]>0)))])
  }
  return(fits)
}


sim_slopes <- list()
for (i in 1:500) {
  slopes <- vector()
  for (j in 1:12) {
    slopes[j] <- sim_fits[[i]][[j]]$coefficients[2]
  }
  sim_slopes[[i]] <- slopes
}

sim_results <- vector()
for (i in 1:500) {
  sim_results[i] <- lm(sim_slopes[[i]][9:12]~neg_stages[90:93])$coefficients[2]
}



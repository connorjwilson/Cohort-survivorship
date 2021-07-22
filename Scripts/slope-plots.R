####################################
# slope-plots.R
####################################
# You must these scripts
# first for this script to work
####################################

source("Scripts/slope-calculation.R")
source("Scripts/sp-slope-calculation.R")

######################################
## Slopes over time comparison plot
######################################

plot(neg_stages[82:93], na_slopes, ylim=c(-1, 0), type="b")
points(neg_stages[82:93], sa_slopes[1:12], col="red", type="b")
points(neg_stages[82:93], eu_slopes, col="blue", type="b")
points(neg_stages[82:93], as_slopes, col="gray", type="b")
points(neg_stages[82:93], af_slopes, col="darkblue", type="b")
points(neg_stages[82:93], au_slopes, col="darkgreen", type="b")

##################################
## Slopes over time (Sp. level)
##################################

plot(neg_stages[82:93], nasp_slopes, ylim=c(-1,0), type="b")
points(neg_stages[82:93], eusp_slopes, col="blue", type="b")
points(neg_stages[82:93], assp_slopes, col="gray", type="b")

####################################
## Species/genus level comparison
####################################

plot(neg_stages[82:93], na_slopes, ylim=c(-1, 0), type="b")
points(neg_stages[82:93], nasp_slopes, col="red", type="b")

plot(neg_stages[82:93], eu_slopes, ylim=c(-1, 0), type="b")
points(neg_stages[82:93], eusp_slopes, col="red", type="b")

plot(neg_stages[82:93], as_slopes, ylim=c(-1, 0), type="b")
points(neg_stages[82:93], assp_slopes, col="red", type="b")


################################################
## Plotting survivorship curves vs simulations
################################################

library(paleotree)
library(divDyn)
data(stages)

load("Data/simulations.RData")
sim_data <- simulations
rm(simulations)

simu_ranges <- list()
for (i in 1:500) {
  simu_ranges[[i]] <- as.data.frame(fossilRecord2fossilRanges(sim_data[[i]], 
                                                              ranges.only=TRUE))
}

simu_ranges_notna <- list()
for (i in 1:500) {
  simu_ranges_notna[[i]] <- simu_ranges[[i]][!is.na(simu_ranges[[i]]$FAD),]
}

classify <- function (input) {
  output <- input
  
  for (i in 1:nrow(input)) {
    output[i,1] <- which((input[i,1]<=stages$bottom & input[i,1]>=stages$top)==TRUE)
    output[i,2] <- which((input[i,2]<=stages$bottom & input[i,2]>=stages$top)==TRUE)
  }
  
  return(output)
}

simu_ranges_stg <- list()
for (i in 1:500) {
  simu_ranges_stg[[i]] <- classify(simu_ranges_notna[[i]])
}

for (i in 1:500) {
  simu_ranges_stg[[i]]$duration <- (simu_ranges_stg[[i]]$LAD - simu_ranges_stg[[i]]$FAD)
  
}

simu_surv <- list()
for (i in 1:500) {
  simu_surv[[i]] <- survivors(fl=simu_ranges_stg[[i]])
}

tsplot(stages, shading="series", boxes="sys",
       xlim=c(66,0), ylab="proportion of survivors present",
       ylim=c(0.01,1),plot.args=list(log="y"))
for (j in 1:150) {for(i in 1:ncol(simu_surv[[j]])) lines(stages$mid, simu_surv[[j]][,i])}

for(i in 1:ncol(surv)) lines(stages$mid, surv[,i], col="blue", lwd=3)
for(i in 1:ncol(eu_surv)) lines(stages$mid, eu_surv[,i], col="red", lwd=3)

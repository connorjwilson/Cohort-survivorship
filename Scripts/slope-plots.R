####################################
# slope-plots.R
####################################
# You must run slope-calculation.R
# first for this script to work
####################################

source("Scripts/slope-calculation.R")


######################################
## Slopes over time comparison plot
######################################

plot(neg_stages[82:93], na_slopes, ylim=c(-1, 0), type="b")
points(neg_stages[82:93], sa_slopes[1:12], col="red", type="b")
points(neg_stages[82:93], eu_slopes, col="blue", type="b")
points(neg_stages[82:93], as_slopes, col="gray", type="b")


####################################
## Species/genus level comparison
####################################

sasp_data <- read.csv("Data/pbdb_sa_sp.csv", skip = 20, header = TRUE)

stgMin <- categorize(sasp_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(sasp_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

sasp_data$stg <- rep(NA, nrow(sasp_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
sasp_data$stg[stgCondition] <- stgMin[stgCondition]

sasp_surv <- survivors(sasp_data, bin='stg', tax='accepted_name')
sasp_cohort <- sasp_surv[82:95,82:95]

sasp_fits <- list()
for (i in 1:14) {
  sasp_fits[[i]] <- lm(log(sasp_cohort[i:max(which(sasp_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(sasp_cohort[,i]>0)))])
}

sasp_slopes <- vector()
for (i in 1:12) {
  sasp_slopes[i] <- sasp_fits[[i]]$coefficients[2]
}

sasp_result <- lm(sasp_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], sasp_slopes, type="b")
abline(sasp_result)

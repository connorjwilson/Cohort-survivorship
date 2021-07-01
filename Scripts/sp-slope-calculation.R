############################
# sp-slope-calculation.R
############################

library(divDyn)
data(keys)
data(stages)

neg_stages <- -1 * stages$mid

########################
## South America
########################

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


######################
## North America
######################

nasp_data <- read.csv("Data/pbdb_na_sp.csv", skip = 20, header = TRUE)

stgMin <- categorize(nasp_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(nasp_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

nasp_data$stg <- rep(NA, nrow(nasp_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
nasp_data$stg[stgCondition] <- stgMin[stgCondition]

nasp_surv <- survivors(nasp_data, bin='stg', tax='accepted_name')
nasp_cohort <- nasp_surv[82:95,82:95]

nasp_fits <- list()
for (i in 1:14) {
  nasp_fits[[i]] <- lm(log(nasp_cohort[i:max(which(nasp_cohort[,i]>0)),i])
                       ~ neg_stages[(81+i):(81+max(which(nasp_cohort[,i]>0)))])
}

nasp_slopes <- vector()
for (i in 1:12) {
  nasp_slopes[i] <- nasp_fits[[i]]$coefficients[2]
}

nasp_result <- lm(nasp_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], nasp_slopes, type="b")
abline(nasp_result)


#######################
## Europe
#######################

eusp_data <- read.csv("Data/pbdb_eu_sp.csv", skip = 20, header = TRUE)

stgMin <- categorize(eusp_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(eusp_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

eusp_data$stg <- rep(NA, nrow(eusp_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
eusp_data$stg[stgCondition] <- stgMin[stgCondition]

eusp_surv <- survivors(eusp_data, bin='stg', tax='accepted_name')
eusp_cohort <- eusp_surv[82:95,82:95]

eusp_fits <- list()
for (i in 1:14) {
  eusp_fits[[i]] <- lm(log(eusp_cohort[i:max(which(eusp_cohort[,i]>0)),i])
                       ~ neg_stages[(81+i):(81+max(which(eusp_cohort[,i]>0)))])
}

eusp_slopes <- vector()
for (i in 1:12) {
  eusp_slopes[i] <- eusp_fits[[i]]$coefficients[2]
}

eusp_result <- lm(eusp_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], eusp_slopes, type="b")
abline(eusp_result)


#######################
## Asia
############
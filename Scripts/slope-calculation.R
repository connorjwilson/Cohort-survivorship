#########################
# slope-calculation.R
#########################

library(divDyn)
data(keys)
data(stages)

neg_stages <- -1 * stages$mid

####################
## North America
####################

na_data <- read.csv("Data/pbdb_na.csv", skip = 20, header = TRUE)

stgMin <- categorize(na_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(na_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

na_data$stg <- rep(NA, nrow(na_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
na_data$stg[stgCondition] <- stgMin[stgCondition]

surv <- survivors(na_data, bin='stg', tax='accepted_name')
na_cohort <- surv[82:95,82:95]

na_fits <- list()
for (i in 1:14) {
  na_fits[[i]] <- lm(log(na_cohort[i:max(which(na_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(na_cohort[,i]>0)))])
}

na_slopes <- vector()
for (i in 1:12) {
  na_slopes[i] <- na_fits[[i]]$coefficients[2]
}


na_result <- lm(na_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], na_slopes)
abline(na_result)

### Plots Slope of each cohort

# pdf("Cohorts.pdf",height = 20, width = 3)
# par(mfrow = c(13,1), mar = c(3, 3, 0, 0.5))
# for (i in 1:13) {
#   plot(log(na_cohort[i:max(which(na_cohort[,i]>0)),i])
#        ~ neg_stages[(81+i):(81+max(which(na_cohort[,i]>0)))], 
#        ylim = c(-6,0), xlim = c(-70,0))
#   abline(na_fits[[i]])
# }
# dev.off()

# plot(stages[82:95,9], na_slopes)
# plot(neg_stages[82:95], stages[82:95,9])

####################
## South America
####################

sa_data <- read.csv("Data/pbdb_sa.csv", skip = 21, header = TRUE)

stgMin <- categorize(sa_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(sa_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

sa_data$stg <- rep(NA, nrow(sa_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
sa_data$stg[stgCondition] <- stgMin[stgCondition]

sa_surv <- survivors(sa_data, bin='stg', tax='accepted_name')
sa_cohort <- sa_surv[82:95,82:95]

sa_fits <- list()
for (i in 1:14) {
  sa_fits[[i]] <- lm(log(sa_cohort[i:max(which(sa_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(sa_cohort[,i]>0)))])
}

sa_slopes <- vector()
for (i in 1:12) {
  sa_slopes[i] <- sa_fits[[i]]$coefficients[2]
}

sa_result <- lm (sa_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], sa_slopes, type="b")
abline(sa_result)

##############
## Europe
##############

eu_data <- read.csv("Data/pbdb_eu.csv", skip = 21, header = TRUE)

stgMin <- categorize(eu_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(eu_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

eu_data$stg <- rep(NA, nrow(eu_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
eu_data$stg[stgCondition] <- stgMin[stgCondition]

eu_surv <- survivors(eu_data, bin='stg', tax='accepted_name')
eu_cohort <- eu_surv[82:95,82:95]

eu_fits <- list()
for (i in 1:14) {
  eu_fits[[i]] <- lm(log(eu_cohort[i:max(which(eu_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(eu_cohort[,i]>0)))])
}

eu_slopes <- vector()
for (i in 1:12) {
  eu_slopes[i] <- eu_fits[[i]]$coefficients[2]
}

eu_result <- lm(eu_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], eu_slopes)
abline(eu_result)

###################
## Africa
###################

af_data <- read.csv("Data/pbdb_af.csv", skip = 21, header = TRUE)

stgMin <- categorize(af_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(af_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

af_data$stg <- rep(NA, nrow(af_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
af_data$stg[stgCondition] <- stgMin[stgCondition]

af_surv <- survivors(af_data, bin='stg', tax='accepted_name')
af_cohort <- af_surv[82:95,82:95]

af_fits <- list()
for (i in 2:14) {
  af_fits[[i]] <- lm(log(af_cohort[i:max(which(af_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(af_cohort[,i]>0)))])
}

af_slopes <- vector()
for (i in 2:12) {
  af_slopes[i] <- af_fits[[i]]$coefficients[2]
}

af_result <- lm(af_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], af_slopes, type="b")
abline(af_result)

###############
## Asia
###############

as_data <- read.csv("Data/pbdb_as.csv", skip = 21, header = TRUE)

stgMin <- categorize(as_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(as_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

as_data$stg <- rep(NA, nrow(as_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
as_data$stg[stgCondition] <- stgMin[stgCondition]

as_surv <- survivors(as_data, bin='stg', tax='accepted_name')
as_cohort <- as_surv[82:95,82:95]

as_fits <- list()
for (i in 1:14) {
  as_fits[[i]] <- lm(log(as_cohort[i:max(which(as_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(as_cohort[,i]>0)))])
}

as_slopes <- vector()
for (i in 1:12) {
  as_slopes[i] <- as_fits[[i]]$coefficients[2]
}

as_result <- lm(as_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], as_slopes)
abline(as_result)

##################
## Australia
##################

au_data <- read.csv("Data/pbdb_au.csv", skip = 21, header = TRUE)

stgMin <- categorize(au_data[,"early_interval"], keys$stgInt)
stgMax <- categorize(au_data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

au_data$stg <- rep(NA, nrow(au_data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
au_data$stg[stgCondition] <- stgMin[stgCondition]

au_surv <- survivors(au_data, bin='stg', tax='accepted_name')
au_cohort <- au_surv[82:95,82:95]

au_fits <- list()
for (i in 1:14) {
  if (i %in% c(1, 2, 4, 5, 6)) {
    au_fits[[i]] <- NA
  } else {
  au_fits[[i]] <- lm(log(au_cohort[i:max(which(au_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(au_cohort[,i]>0)))])
  }
}

au_slopes <- vector()
for (i in 1:12) {
  if (!is.na(au_fits[[i]]))
  au_slopes[i] <- au_fits[[i]]$coefficients[2]
}

au_result <- lm(au_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], au_slopes)
abline(au_result)

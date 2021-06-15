library(divDyn)
data(keys)
data(stages)

data <- read.csv("Data/pbdb_na.csv", skip = 20, header = TRUE)

neg_stages <- -1 * stages$mid

## North America
na_cohort <- surv[82:95,82:95]

na_fits <- list()
for (i in 1:14) {
  na_fits[[i]] <- lm(log(na_cohort[i:max(which(na_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(na_cohort[,i]>0)))])
}

na_slopes <- vector()
for (i in 1:14) {
  na_slopes[i] <- na_fits[[i]]$coefficients[1]
}

na_result <- lm (na_slopes ~ neg_stages[82:95])
plot(neg_stages[82:95], na_slopes)
abline(na_result)

plot(stages[82:95,9], na_slopes)
plot(neg_stages[82:95], stages[82:95,9])

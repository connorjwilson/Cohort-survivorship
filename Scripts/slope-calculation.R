neg_stages <- -1 * stages$mid

## North America
na_cohort <- surv[82:95,82:95]

na_fit1 <- lm(na_cohort[,1] ~ stages$mid[82:95])

na_fit2 <- lm(log(na_cohort[1:6,1]) ~ neg_stages[82:87])
na_slopes <- vector()

neg_stages <- -1 * stages$mid

## North America
na_cohort <- surv[82:95,82:95]

na_slopes <- vector()
for (i in 1:14){
  na_slopes[i] <- lm(log(na_cohort[i:max(which(na_cohort[,i]>0)),i])
                     ~ neg_stages[81+i:81+max(which(na_cohort[,i]>0))])
}

na_fit1 <- lm(na_cohort[,1] ~ stages$mid[81+i:81+max(which(na_cohort[,i]>0))])

na_fit2 <- lm(log(na_cohort[1:max(which(na_cohort[,i]>0)),1]) ~ neg_stages[82:87])


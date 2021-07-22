library(paleotree)

test <- simFossilRecord(p = 0.25, q = 0.25, r = 0.1, startTaxa = 100, 
                        totalTime = 66, nTotalTaxa = c(0,Inf), 
                        nExtant = c(0,Inf), nSamp = c(0,Inf), print.runs = T, 
                        plot = T)

test <- simFossilRecord(p = 0.25, q = 0.25, r = "0.902 + 0.006059 * T + 0.00006036 * T^2",
                        startTaxa = 100, 
                        totalTime = 66, nTotalTaxa = c(0,Inf), 
                        nExtant = c(0,Inf), nSamp = c(0,Inf), print.runs = T, 
                        plot = T)


ranges_data <- fossilRecord2fossilRanges(test, ranges.only=TRUE)


########################
library(dplyr)
library(tidyr)


test_fl <- fadlad(data, 'accepted_name', 'stg')
########################

ranges_df <- as.data.frame(ranges_data)
ranges_df1 <- ranges_df[!is.na(ranges_df$FAD),]

# classify <- function (input){
# 
# output <- input
#       
#   for (i in 1:nrow(input)) {
#       for (j in 1:nrow(stages)) {
#       
#           if (input[i,1]<=stages$bottom[j] & input[i,1]>stages$top) {
#               output[i,1] <- stages$stg[j]
#           }
#         
#           if (input[i,2]<=stages$bottom[j] & input[i,2]>stages$top) {
#               output[i,2] <- stages$stg[j]
#           }
#     }
#   }
#   return(output)
# }

ranges_cf <- classify(ranges_df1)

classify <- function (input){
  
  output <- input
  
  for (i in 1:nrow(input)) {
    output[i,1] <- which((input[i,1]<=stages$bottom & input[i,1]>=stages$top)==TRUE)
    output[i,2] <- which((input[i,2]<=stages$bottom & input[i,2]>=stages$top)==TRUE)
  }
  return(output)
}

ranges_cf$duration <- ranges_cf$LAD - ranges_cf$FAD

sim_surv <- survivors(fl=ranges_cf)
tsplot(stages, shading="series", boxes="sys",
       xlim=c(66,0), ylab="proportion of survivors present",
       ylim=c(0.01,1),plot.args=list(log="y"))
for(i in 1:ncol(sim_surv)) lines(stages$mid, sim_surv[,i])

sim_cohort <- sim_surv[82:95, 82:95]

sim_fits <- list()
for (i in 1:14) {
  sim_fits[[i]] <- lm(log(sim_cohort[i:max(which(sim_cohort[,i]>0)),i])
                     ~ neg_stages[(81+i):(81+max(which(sim_cohort[,i]>0)))])
}

sim_slopes <- vector()
for (i in 1:12) {
  sim_slopes[i] <- sim_fits[[i]]$coefficients[2]
}


sim_result <- lm(sim_slopes ~ neg_stages[82:93])
plot(neg_stages[82:93], sim_slopes, type="b")
abline(sim_result)

library(paleotree)

test <- simFossilRecord(p = 0.25, q = 0.25, r = 0.1, startTaxa = 100, 
                        totalTime = 66, nTotalTaxa = c(0,Inf), 
                        nExtant = c(0,Inf), nSamp = c(0,Inf), print.runs = T, 
                        plot = T)

ranges_data <- fossilRecord2fossilRanges(test, ranges.only=TRUE)


########################
library(dplyr)
library(tidyr)


test_fl <- fadlad(data, 'accepted_name', 'stg')
########################


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


library(divDyn)
data(keys)
data(stages)

data <- read.csv("Data/pbdb_na.csv", skip = 20, header = TRUE)

stgMin <- categorize(data[,"early_interval"], keys$stgInt)
stgMax <- categorize(data[,"late_interval"], keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

data$stg <- rep(NA, nrow(data))
stgCondition <- c(which(stgMax==stgMin), which(stgMax==-1))
data$stg[stgCondition] <- stgMin[stgCondition]

results <- divDyn(data, bin="stg", tax="accepted_name")

tsplot(stages, boxes="sys", shading="sys", xlim=82:95, ylim=c(0, 1000),
       ylab="richness", xlab="age (Ma)")
lines(stages$mid, results$divCSIB, col="blue", lwd=2)


## Occurences
sampStg <- binstat(data, tax="accepted_name", bin="stg", coll="collection_no", duplicates=FALSE)

tsplot(stages, boxes="sys", shading="sys", xlim=83:95, ylim=c(0,5000),
       ylab="Number of entries", xlab="Age (Ma)")
lines(stages$mid, sampStg$occs, lwd=2)
lines(stages$mid, sampStg$colls, lwd=2, col="blue")
legend("top", bg="white", legend=c("occurrences", "collections"),
       col=c("black", "blue"), lwd=2, inset=c(0.15,0.01), cex=1)

## Survivorship
surv <- survivors(data, bin='stg', tax='accepted_name')
tsplot(stages, shading="series", boxes="sys",
       xlim=c(66,0), ylab="proportion of survivors present",
       ylim=c(0.01,1),plot.args=list(log="y"))
for(i in 1:ncol(surv)) lines(stages$mid, surv[,i])

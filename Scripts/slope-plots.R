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
points(neg_stages[82:93], af_slopes, col="darkblue", type="b")


####################################
## Species/genus level comparison
####################################



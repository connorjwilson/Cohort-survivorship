########################
# parameter-fitting.R
########################

slopes <- c(na_slopes)
times <- c(neg_stages[82:93])

df <- as.data.frame(cbind(slopes, times))

plot(times, slopes)
fitting <- lm(slopes ~ times+ I(times^2))

x <- with(df, seq(min(times), max(times), length.out=2000))
y <- predict(fitting, newdata=data.frame(times=x))

lines(x,y)

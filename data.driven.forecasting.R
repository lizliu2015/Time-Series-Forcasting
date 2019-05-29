setwd("/Users/Liz/...")
rm(list=ls(all=TRUE)) #clear data

data<-read.csv("data_week.csv", sep=",",dec=".",header=T) 	# weekly data
names(data)
head(data,5)

sales = data[,3]/1000000		# millions in xyz currency
temp = data[,11]				# centigrade



##### Representing Data as Time Series Objects #####

yy = ts(sales, frequency = 52,start = c(2015,1)) # coverts sales data as time series object with start date and weekly frequency
plot.ts(yy)									# ALWAYS plot time series to see patterns: trend, cycle, variance over time

xx = ts(temp, frequency = 52, start = c(2015,1))
plot.ts(xx)



##### Moving Average #####

library("TTR")		# install TTR package first and then load the library
yy4 = SMA(yy,4)			# 4-wekk moving average
yy13 = SMA(yy,13)		# 13 week (quarterly) moving average
yy52 = SMA(yy,52)		# annual moving average

plot.ts(cbind(yy,yy4,yy13,yy52))		# understand the "smoothing" concept due to the averaging operation


##### Time Series Decomposition  #####

## Sales: What's the growth rate? How do we obtain seasonally adjusted sales?
yd = decompose(yy) 
yd.trend = yd$trend
yd.seasonal = yd$seasonal
yd.random = yd$random
yy.season.adj = yy - yd.seasonal									# seasonally adjusted sales
plot.ts(cbind(yy,yd.trend, yd.seasonal, yd.random, yy.season.adj))

## Temperature: Are temperatures increasing?
xd = decompose(xx) 
xd.trend = xd$trend
xd.seasonal = xd$seasonal
xd.random = xd$random
plot.ts(cbind(xx,xd.trend, xd.seasonal, xd.random))



##### Simple Exponential Smoothing #####

out1 = HoltWinters(xx, beta=TRUE, gamma=TRUE) 				# Holt-Winters Filtering - only level updating b/c beta and gamma are zeros
out1																# output -- see alpha estimate
out1$fitted 														# fitted values in training data
plot(out1)															# graph of actual (black) vs fitted (red)

##### Out of Sample Forecasts

library("forecast")								# install "forecast" package and load the library
out2 = forecast:::forecast.HoltWinters(out1, h = 26, level = c(68, 95))	
# forecast horizon 26 weeks. CI 50%, 80%, 95% 

# something wierd as should directly call forecast.HoltWinters() function, but does not. Hence, to force it, I used the prefix forecast::: 
plot(out2)


##### Autocorrelation Function (ACF)

## Can forecasts be improved upon using information in the residuals?
## To answer this Q, we need to assess whether the residuals contain information. That is, are residuals correlated at various time lags?

resid = out2$residuals				# output the residuals
corr.lag = Acf(resid[2:130])		# all correlations insignificant ==> no information in residuals ==> Model cannot be improved further

# Remarks
# First obs is NA and hence eliminated.
# Acf() gives autocorrelation function w/o the correlation = 1 at lag = 0. 
# acf() with the lower case "a" gives correlations at all lags, including lag 0, which by definition is 1. 
# acf()tends to mask small correlations in the plot due to the wider y-scale to include 1.
# The confidence interval at 95% is 2/sqrt(T). Here T = 130 ==> 0.175


##### Partial Autocorrelation Function (PACF)

corr.lag = Pacf(resid[2:130]) 		# "Partial" refers to correlations obtained *after controlling for the effects of previous correlations*


# Remarks
# Plotting and then assessing whether correlations spike above confidence bands requires manual process. 
# To automate this assessment, we use Ljung-Box test. If p-value < 0.05, then significant correlations exist; else, not



#### Ljung-Box test 
Box.test(resid, type = "Ljung-Box")		

# The above test yeilds a large p-value ==> so no significant correlations as we verified visually from acf and pacf






  
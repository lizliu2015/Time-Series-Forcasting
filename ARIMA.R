setwd("/Users/Liz/.../")
rm(list=ls(all=TRUE)) 	# clear data
install.packages('forecast')
library("forecast")
install.packages('tseries')
library("tseries") 		# reqired for adf.test of stationarity

data<-read.csv("data_week.csv", sep=",",dec=".",header=T) 	# weekly data
names(data)
head(data,5)

# getting 3 variables - avg_cloud_index, avg_temp, and avg_hours_sun from dataset
cloud = data[,10]	    #index
temp = data[,11]			#degree
sun = data[,14]       #hour

X = cbind(cloud, temp, sun)
X
# PCA dimensionality reduction
library(stats)
X.pca <- prcomp(X, scale. = TRUE, rank.=1, retx=True)
summary(X.pca)
X.pca
pc <- X.pca$x

##### Representing Data as Time Series Object #####

yy = ts(pc, frequency = 52,start = c(2015,1))		# coverts data as time series object with start date and frequency (weekly here)
plot.ts(yy)									              			# ALWAYS plot time series to see patterns: we can observe seasonality in the plot

##### General Process for Fitting ARIMA(p,d,q) x (P, D, Q) Models #####

# Typical values of (p,q) are 0, 1, 2. So you will have 3 ps x 3 qs = 9 models

# ARMA models require STATIONARY time series as inputs. Stationarity implies mean and variance are approx constant over time)

# To assess stationarity, use adf.test(). If series is non-stationary, then take first difference. If first-differenced series still not stationary, take higher order difference. Usually d = 2 suffices for most time series. 

# Given p, d, q taking 3 values (0,1,2), you will have a set of 27 models. Apply AICC to select the best one or a set of few good ones. 

# If seasonal arima with (P, D, Q) were also included, then you will have additonal 27 models. So a total of 54 models.  Apply AICC to select the best one or a set of few good ones. 
 


## Let's learn the process using sales series

## Step 1. Is the time series stationary? 

# Use Augmented Dickey-Fuller Test to test stationarity == > large p-value means nonstationary

plot.ts(yy)								# ALWAYS plot series to eyeball the presence of trends and/or variability over disjoint regimes

# install and load "tseries" package 
adf.test(yy)							# formal test for stationarity ==> if p-value is large (> 0.10), then nonstationary
                          ## p-value = 0.2487 > 0.10, nonstationary

yd = diff(yy,differences = 1)			
plot.ts(yd)								# looks stationary visually
adf.test(yd)							# estimated p = 0.01 ==> small pvalue (< 0.10) ==> so yd is stationary ==> fix d = 1 in ARIMA models to be fitted
                          ## p-value < 0.01,  first order differencing is stationary



## Step 2. Decide AR(p) or MA(q) or both ARMA(p,q). Use the stationary series from Step 1. 

# To decide AR (p), plot Pacf. AR(p) signature = Pacf becomes zero at some lag p

pacf(yd, lag.max = 10)					# Pacf suggests p = 3 


# To decide MA, plot Acf. MA(q) signature = Acf becomes zero at some lag q

acf(yd, lag.max = 10)				# # Acf suggests q = 2 




## Step 3. Fit several ARIMA models and calculate the information criteria
aic_score = c()
aicc_score = c()
bic_score = c()
p = ncol(yy)
n = nrow(yy)

# p ranges from 0-3, d = 1, q ranges from 0-2

for (i in 0:3) {
  for (j in 0:2) {
    model = arima(yy, order=c(i,1,j))			
    aic = model$aic
    aicc = aic + 2*(p+1)*(p+2)/(n-p-2)
    bic = aic - 2*p + n*log(p)
    aic_score <- c(aic_score, aic)
    aicc_score <- c(aicc_score, aicc)
    bic_score <- c(bic_score, bic)
  }
}
  
scores = cbind(aic_score, aicc_score, bic_score)
scores  

# AIC --- similar to R^2 
#     --- 5% < p/n < 15%
# AICC = AIC + 2(p+1)(p+2)/(n-p-2) ???--- similar to adjusted R^2, penalize for the complexity of model
#                                  --- if small sample, large p, p/n > 15%, use AICc
# BIC = n*Ln(sigma^2) + n*Ln(p)    --- good when n is much larger than p, p/n < 5% BIC


# model 2, 10, 12 are the best


# auto.arima 
model_auto= auto.arima(yy)		# fits ARIMA(p,d,q) x (P, D, Q) automatically
summary(model_auto)

# ARIMA(0,0,0)(0,1,0)[52] 
# AIC=267.68   AICc=267.74   BIC=270.04


model_auto_1= auto.arima(yy, seasonal = FALSE)		# if restrict to non-seasonal model
summary(model_auto_1)
#ARIMA(1,0,1) with zero mean 


# Consider Seasonal ARIMA(p,d,q) x (P, D, Q) components when seasonality is expected/suspected
aic_score = c()
aicc_score = c()
bic_score = c()
p = ncol(yy)
n = nrow(yy)

      for (a in 0:1) {
        for (b in 0:1) {
          for (c in 0:1) {
            model = arima(yy, order=c(1, 0, 1), seasonal = list(order = c(a,b,c), period = 52))			
            aic = model$aic
            aicc = aic + 2*(p+1)*(p+2)/(n-p-2)
            bic = aic - 2*p + n*log(p)
            aic_score <- c(aic_score, aic)
            aicc_score <- c(aicc_score, aicc)
            bic_score <- c(bic_score, bic)
          }
        }
      }

scores = cbind(aic_score, aicc_score, bic_score)
scores  
# model 3 and model 7 are the best
# model 3: P, Q, R = 0, 1, 0
# model 7: P, Q, R = 1, 1, 0

# Consider Seasonal ARIMA(p,d,q) x (P, D, Q) components when seasonality is expected/suspected
aic_score = c()
aicc_score = c()
bic_score = c()
p = ncol(yy)
n = nrow(yy)

      for (a in 0:1) {
        for (b in 0:1) {
          for (c in 0:1) {
            model = arima(yy, order=c(0,0,0), seasonal = list(order = c(a,b,c), period = 52))			
            aic = model$aic
            aicc = aic + 2*(p+1)*(p+2)/(n-p-2)
            bic = aic - 2*p + n*log(p)
            aic_score <- c(aic_score, aic)
            aicc_score <- c(aicc_score, aicc)
            bic_score <- c(bic_score, bic)
          }
        }
      }

scores = cbind(aic_score, aicc_score, bic_score)
scores

## model 3, 4, 7 are the best among 8 models
## model 3: P, Q, R = 0, 1, 0
## model 4: P, Q, R = 0, 1, 1
## model 7: P, Q, R = 1, 1, 0

## Step 5. I dentify a **set of few good models**. Good models are within plus/minus 1 point difference from the smallest value of an info criterion.  
m1 = arima(yy, order=c(1,0,1), seasonal = list(order = c(0,1,0), period = 52))	
m2 = arima(yy, order=c(1,0,1), seasonal = list(order = c(0,1,1), period = 52))

m3 = arima(yy, order=c(0,0,0), seasonal = list(order = c(0,1,0), period = 52))	# this one is the best among five
m4 = arima(yy, order=c(0,0,0), seasonal = list(order = c(0,1,1), period = 52))
m5 = arima(yy, order=c(0,0,0), seasonal = list(order = c(1,1,0), period = 52))

## Step 6. Make Out-of-Sample Forecasts with Prediction Interval based on your retained model

m1.predict = forecast:::forecast.Arima(m1, h = 52, level = c(68, 80, 90, 95))
plot(m1.predict)

m2.predict = forecast:::forecast.Arima(m2, h = 52, level = c(68, 80, 90, 95)) 
plot(m2.predict)

m3.predict = forecast:::forecast.Arima(m3, h = 52, level = c(68, 80, 90, 95)) #same as auto.arima model
plot(m3.predict)

m4.predict = forecast:::forecast.Arima(m4, h = 52, level = c(68, 80, 90, 95))
plot(m4.predict)

m5.predict = forecast:::forecast.Arima(m5, h = 52, level = c(68, 80, 90, 95))
plot(m5.predict)
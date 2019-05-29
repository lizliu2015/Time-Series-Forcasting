setwd("C:/Users/Liz/...")
rm(list=ls(all=TRUE)) #clear data
library(ggplot2)
library(gridExtra)
library(numDeriv)
library(Hmisc)
library(GenSA)

##### Read Data #####

data<-read.csv("chocolate.csv", sep=",",dec=".",header=T) 	# weekly data
names(data)
head(data,5)


##### Plot Data #####

p1 = ggplot(data, aes(x=Week, y = GRP)) + geom_bar(stat="identity") 
p2 = ggplot(data, aes(x=Week, y = Awareness)) + geom_line() + geom_point(colour="blue",size=2, shape = 21,fill="white") 

grid.arrange(p1,p2, nrow=2)  # present both GRP and Awareness plots stacked in 2 rows
theme_set(theme_bw())
theme(panel.grid.minor = element_blank()) # turn off minor grid in each plot


##### Variables #####

grp = data[,2]/10			# GRP = reach x frequency = percentage market reached x number of times per week = ad spending intensity
aware = data[,3]			# percentage market aware of ad campaign 
on = ifelse(grp>0,1,0)		# indicator for advertising is on (= 1) or off (= 0)
x1 = data[,4]				# Ad campaign A is on
x2 = data[,5]				# Ad campaign B is on
x3 = data[,6]				# Ad campaign C is on

aware.lag = Lag(aware)
out = lm(aware~aware.lag+grp)

# I ran lm to check the grp effect on awareness -- beta value was 0.056. 
# So I divided awareness by 10 to yield grp effect = 0.5, which is on the same scale as the beta value of lag.aware effect.
# Keeping parameters to be estimated on the same magnitude is a good practice when using nonlinear optimization to get convergence.


ra = nrow(data)							# number of observations


##### KALMAN FILTER LIKELIHOOD #####

neg_likeli = function(b)
{
	# Setup non-time varying system matrices {Z, T, c, d, Q, R} outside the for-loop 

	z = as.matrix(cbind(1,0));  		# Link matrix Z 
 	tt = diag(2)			  			# Transition matrix T
 	ct = 0								# Observation drift vector 
 	dt = matrix(0,2,1)					# Transition drift vector
 	
 	qq = b[1]^2							# Observation noise (scalar b/c of univariate awareness time series)
 	
 	rr = matrix(0,2,2)					# Transition noise
 	rr[1,1] = b[2]^2			
 	rr[2,2] = b[3]^2
 	
   	a0 = b[4:5]							# mean initial state to be estimated
   	p0 = 25*diag(2)      	  			# covariance of initial state (difused prior with var 25 = mean(aware)
   
   	# model parameters (non time varying outside the for-loop)
    delta = b[6]						# forgetting rate
    c = b[7]							# copy wearout
    w = 0								# repetition wearout
    
    # initialize
   	at = a0				  	  
   	pt = p0				  	   
   	likeli= 0				  
   	
   
   	for (ii in 1:ra)		
       { 
       	
       	# Setup time varying system matrices {Z, T, c, d, Q, R} inside this for-loop
       	
       	tt[1,1] = exp(-delta)
       	tt[1,2] = grp[ii]
       	tt[2,1] = 0
       	tt[2,2] = exp(-c-w*grp[ii] - delta*(1-on[ii]))
       	
       	dt[1] = 0
       	dt[2] = delta*(1-on[ii])
       	
       	
       	# Kalman Filter Recursions
       	# Notation ==> xtm1 means x at (t-1). x at t|t-1 denoted as xt.tm1. Read "m"  as minus
       	
	  	atm1 = at			# current period's prior mean is last period's posterior mean		
       	ptm1 = pt 			# current period's prior covar is last period's posterior covar			
       	
       	
       	# TIME UPDATE (based on model dynamics)
       	at.tm1 = tt %*% atm1 + dt				
       	pt.tm1 = tt %*% ptm1 %*% t(tt) + rr  	
       		
       		
       	yhat = z %*% at.tm1 + ct					# predicted awareness
       	err = aware[ii] - yhat						# prediction error 
       	f = z %*% pt.tm1 %*% t(z) + qq				# variance in data 
        kgain = pt.tm1 %*% t(z) %*% solve(f) 		# Kalman Gain Factor
        
        
        # INFORMATION UPDATE (based on observed data)
      	at = at.tm1 + kgain %*% err;				# posterior mean
       	pt = pt.tm1 - kgain %*% z %*% pt.tm1		# posterior variance
      
      
      	# minus log-likelihood 
       	likeli = likeli + 0.5*(log(f) + err^2/f)  
       	
       	ii=ii+1										# process next observation from ii = 1,...,ra
       	}
return(likeli) 										# returns minus log-likelhood value for minimization
}


##### Obtain Good Starting Values Using Derivative-Free Solver like Simulated Annealing) #####

b0 = matrix(0.5,7,1)
sa.out = optim(b0,neg_likeli, method = "SANN", control = list(trace = TRUE, maxit = 5000))



##### MLE: Maximum Likelihood Estimation using BFGS) #####

b1 = sa.out$par  	
result = optim(b1,neg_likeli, method = "BFGS", hessian = TRUE, control = list(trace = TRUE, maxit = 10000))

converge = result$convergence
maxLik = -result$value
est = result$par
se = sqrt(diag(solve(result$hessian)))
tvals = est/se

# MLE Results
out1 = cbind(abs(est),abs(se),tvals)
colnames(out1) <- c("Estimates","Std Errors", "t-values")
rownames(out1) <- c("Awareness Noise", "Goodwill Noise", "Ad Effectiveness Noise", "Initial Value A0", "Initial Value b0", 
"Forgetting Rate, Delta", "Copy Wearout, c")





#####  Means and Confidence Intervals of the Estimated State Vector  #####

state = matrix(0,ra,2)							
lo.state = matrix(0,ra,2)
hi.state = matrix(0,ra,2)



alpha_hat = function(b)
{

	z = as.matrix(cbind(1,0));  			# Link matrix Z 
 	tt = diag(2)			  			# Transition matrix T
 	ct = 0								# Observation drift vector 
 	dt = matrix(0,2,1)					# Transition drift vector
 	
 	qq = b[1]^2							# Observation noise (scalar b/c of univariate awareness time series)
 	
 	rr = matrix(0,2,2)					# Transition noise
 	rr[1,1] = b[2]^2			
 	rr[2,2] = b[3]^2
 	
   	a0 = b[4:5]							# mean initial state to be estimated
   	p0 = 25*diag(2)      	  			# covariance of initial state (difused prior with var 25 = mean(aware)
   
   	
    delta = b[6]						# forgetting rate
    c = b[7]							# copy wearout
    w = 0								# repetition wearout
    
    # initialize
   	at = a0				  	  
   	pt = p0				  	   
   	likeli= 0				  
   	
   
   	for (ii in 1:ra)		
       { 
       	
       	tt[1,1] = exp(-delta)
       	tt[1,2] = grp[ii]
       	tt[2,1] = 0
       	tt[2,2] = exp(-c-w*grp[ii] - delta*(1-on[ii]))
       	
       	dt[1] = 0
       	dt[2] = delta*(1-on[ii])
       	
       	
	  	atm1 = at			# current period's prior mean is last period's posterior mean		
       	ptm1 = pt 			# current period's prior covar is last period's posterior covar			
       	
       	
       	# TIME UPDATE (based on model dynamics)
       	at.tm1 = tt %*% atm1 + dt				
       	pt.tm1 = tt %*% ptm1 %*% t(tt) + rr  	
       		
       		
       	yhat = z %*% at.tm1 + ct						# predicted awareness
       	err = aware[ii] - yhat						# prediction error 
       	f = z %*% pt.tm1 %*% t(z) + qq				# variance in data 
        kgain = pt.tm1 %*% t(z) %*% solve(f) 		# Kalman Gain Factor
        
        
        # INFORMATION UPDATE (based on observed data)
      	at = at.tm1 + kgain %*% err;				# posterior mean
       	pt = pt.tm1 - kgain %*% z %*% pt.tm1		# posterior variance
      
       	# save alpha_hat and lo/hi confidence intervals
       	state[ii,] = t(at) 								# store alpha_hat vector at period ii
      	lo.state[ii,] = t(at - 1.96*sqrt(diag(pt)))		# lo CI 95%
      	hi.state[ii,] = t(at + 1.96*sqrt(diag(pt)))		# hi CI 95%
      	
       	ii=ii+1										# process next observation from ii = 1,...,ra
       	}
	
       	res = cbind(state,lo.state,hi.state)
       	
return(res) 			
}


out2 = cbind(data[,1], alpha_hat(result$par))
colnames(out2) <- c("Week", "Awareness, At","Ad Effectiveness, bt", "Lo Awareness","Lo Beta", "High Awareness","High Beta")





## Plotting Awareness and 95% Confidence Intervals  ## 

df1 <- data.frame(x = 1:ra, ff = out2[,2], L = out2[,4], U = out2[,6])
plot(df1$x, df1$ff, ylim = c(min(cbind(df1$L,0)),max(df1$U)), type = "l", xlab = "Weeks", ylab = "Awareness", main = "Awareness Evolution")

# draws the grey confidence region
polygon(c(df1$x,rev(df1$x)),c(df1$L,rev(df1$U)),col = "grey85", border = FALSE) 
lines(df1$x, df1$ff, lwd = 2)

#add red lines on borders of polygon
lines(df1$x, df1$U, col="red",lty=2)
lines(df1$x, df1$L, col="red",lty=2)

par(new = T)										# suppreses overwriting of the previous plot
points(aware, col = "black", pch = "o")				# overlays points on the previous plot
par(new = F)


## Plotting Ad Wearout and Restoration and 95% Confidence Intervals  ## 

df1 <- data.frame(x = 3:ra, ff = out2[3:ra,3], L = out2[3:ra,5], U = out2[3:ra,7])
plot(df1$x, df1$ff, ylim = c(min(cbind(df1$L,0)),max(df1$U)), type = "l", xlab = "Weeks", ylab = "GRP Effect", main = "Ad Wearout & Restoration")

# draws the grey confidence region
polygon(c(df1$x,rev(df1$x)),c(df1$L,rev(df1$U)),col = "grey85", border = FALSE) 
lines(df1$x, df1$ff, lwd = 2)

#add red lines on borders of polygon
lines(df1$x, df1$U, col="red",lty=2)
lines(df1$x, df1$L, col="red",lty=2)







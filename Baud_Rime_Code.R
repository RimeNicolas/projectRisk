install.packages("evd");install.packages("evdbayes");install.packages("coda");install.packages("ismev")
library(evd);library(evdbayes);library(coda);library(ismev)
#install.packages("reliaR");library(reliaR)
#install.packages("extRemes");library(extRemes)

# Install packages to simplify life
install.packages("hydroTSM");
library(hydroTSM)

# Load data and transform it 
cape = get(load("./data/CAPE_Baud_Rime.RData"))
srh = get(load("./data/SRH_Baud_Rime.RData"))
enso = get(load("./data/NINO34.RData"))
prod = sqrt(cape)*srh

# Define variables
ppd <- 8 #point per day
years <- 37
data_size <- ppd*365*years

# Max of 8 points per day to 1 point per day and add date sequence 
td <- seq(as.Date("1979/01/01"), as.Date("2015/12/31"), "days")
td <- td[!grepl(x = td, pattern = "-02-29$")]
prod_day <- apply(matrix(prod,ncol=ppd,byrow=TRUE),1,max)
prod_day <- zoo(x = prod_day, order.by = td)
# Daily to monthly data taking max per month
prod_m <- daily2monthly(prod_day, FUN=max, na.rm=TRUE) #prod max per month

# Plot daily maximum
plot(td, prod_day)

# Convert zoo to numeric
x_day <- as.numeric(prod_day)
x_m <- as.numeric(prod_m)

# Plot the monthly maximum
plot(x_m)
plot(enso)
plot(x_m, enso)

# QQplot to get an estimate of the parameters
months <- length(x_m)
qqplot(qgumbel(c(1:months)/(months+1)),x_m)
qqline(x_m,distribution=qgumbel)

# Fit GEV 
fit_gev = fgev(x_m)
estimators = fitted(fit_gev) #estimation of the parameters
std = std.errors(fit_gev) #std of the parameters (normal dist.)
par(mfrow=c(2,2))
#par(mfrow=c(1,1))

# Plot pp plot, qq plot, density plot and return level
#pp plot: x: F_GEV(data_i) y: i/(n+1) better for the bulk of the distribution
#qq plot: x: F_GEV^(-1)(i/(n+1)) y: data_i  better for the tail of the distribution
#density plot: x: data_i y: f_GEV^(-1)(i/(n+1))
#return level: x: T y: y_T such that 1 - GEV(y_T) = 1/T 
plot(fit_gev) 

# Plot profile log likelihood for a better (asymmetric) std approximate
plot(profile(fit_gev))

#######################################################
#MCMC
init <- c(5000,1000,0)
mat = diag(c(1e11,1e8,1e8))
pn = prior.norm(mean=c(0,0,0),cov=mat) # normal distribution

psd = c(800,0.10,0.12) # sd for proposal dist (normal dist)
# initial value, prior distribution: normal, likelihood: gev, psd: 
post = posterior(5000,init=init,prior=pn,lh="gev",data=x_m,psd=psd)

#mc = mcmc(post[-c(1:300),])
mc = mcmc(post)
plot(mc)
attr(mc,"ar")

MCMC2<-mcmc(post[-c(1:300),])
acf(MCMC2)

MCMC3<-mcmc(post[500 + 7.5*c(1:337),])
acf(MCMC3)
apply(MCMC3,2,mean)
apply(MCMC3,2,sd)

# Return level MCMC
hist(as.numeric(MCMC3[,3]),nclass=20,prob=T,main="Histogram of xi",xlab="xi")
u.10<-mc.quant(MCMC3,p=0.9,lh="gev")
u.100<-mc.quant(MCMC3,p=0.99,lh="gev")
hist(u.10,nclass=20,prob=T,xlab="10-year return level")
hist(u.100,nclass=20,prob=T,xlab="100-year return level")

###############################################3
# r largest statistics per month
# r = 1, meaning one per month, in order to check with gev and mcmc done previously
months_ <- c(1:months)
x_m <- matrix(x_m,ncol=1)
fit0<-rlarg.fit(x_m) # constant
fit1<-rlarg.fit(x_m,ydat=matrix(months_,ncol=1),mul=c(1)) # linear

# r largest statistics per year (point 5, otherwise only 37 points, not enough)
# r = 1
time_ <- c(1:years)
x1 <- apply(matrix(prod,ncol=(ppd*365),byrow=TRUE),1,max)
x1 <- t(matrix(x1,ncol=(years)))
x1_check <- apply(matrix(x_m,nrow=37,byrow=TRUE),1,max)#check x1 with previous variable x_m

plot(time_,x1)
plot(time_,x1_check)

fit0<-rlarg.fit(x1) # constant
fit1<-rlarg.fit(x1,ydat=matrix(time_,ncol=1),mul=c(1)) # linear

# r = 2
r_ = 2
time2_ <- rep(1:years,each=r_)
x_ <- matrix(prod,ncol=(ppd*365),byrow=TRUE)
x2 <- matrix(nrow = years, ncol = r_)
for (i in c(1:37)){
  x2[i,] = rev(tail(sort(x_[i,], decreasing=FALSE),r_)) 
  #print(x2[i,])
}
plot(time2_,t(x2))

fit0<-rlarg.fit(x2) # constant
fit1<-rlarg.fit(x2,ydat=matrix(time_,ncol=r_),mul=c(1)) # linear

# r = 12 same number of points than first part
r_ = 12
time12_ <- rep(1:years,each=r_)
x_ <- matrix(prod,ncol=(ppd*365),byrow=TRUE)
x12 <- matrix(nrow = years, ncol = r_)
for (i in c(1:37)){
  x12[i,] = rev(tail(sort(x_[i,]),r_)) 
  print(x12[i,])
}

plot(time12_,t(x12))
fit0<-rlarg.fit(x12) # constant
fit1<-rlarg.fit(x12,ydat=matrix(time_,ncol=r_),mul=c(1)) # linear

# r= 50
r_ = 50
time_ <- rep(1:years,each=r_)
x_ <- matrix(prod,ncol=(ppd*365),byrow=TRUE)
x50 <- matrix(nrow = years, ncol = r_)
for (i in c(1:37)){
  x50[i,] = rev(tail(sort(x_[i,]),r_)) 
  #print(x50[i,])
}

plot(time_,t(x50))

# Look if the data is stationary
par(mfrow=c(2,2))
plot(x_m)
plot(time_,x1)
plot(time2_,t(x2))
plot(time12_,t(x12))
par(mfrow=c(1,1))
# seems stationary

#compute correlation prod enso
mat1 <- cbind(x_m,enso)
plot(mat1)
cor(mat1)
# enso in absolute norm
mat2 <- cbind(x_m,abs(enso))
plot(mat2)
cor(mat2)

########################################################################
# Peaks over threshold
plot(prod)
qu.min <- quantile(prod, 0.5) # median value
qu.max <- quantile(prod,(length(prod)-30)/length(prod))
print(paste0("median: ",qu.min,", quantile max: ",qu.max))
mrlplot(prod, tlim=c(qu.min, qu.max))
par(mfrow=c(1,2))
tcplot(prod,tlim=c(qu.min, qu.max))
#or use ismev equivalent
#gpd.fitrange(prod,umin=qu.min, umax=qu.max)

# Choose threshold by hand giving a good variance bias trade-off
th <- 1.8e4
points_year <- 8*365
fit<-fpot(prod,threshold=th,npp=points_year)
par(mfrow=c(2,2))
plot(fit)
# 1.8e4 seems to be the a good threshold, 1.9e4 is less precise
th2 <- 1.9e4
fit2<-fpot(prod,threshold=th2,npp=points_year)
par(mfrow=c(2,2))
plot(fit2)

estimatorsPOT = fitted(fit) #estimation of the parameters
stdPOT = std.errors(fit) #std of the parameters (normal dist.)

# Profile log likelihood 
par(mfrow=c(1,2))
plot(profile(fit))
abline(v=0,col=2,lty=2)
# Zero shape should be considered
fit.gum<-fpot(prod, threshold=th, npp=points_year, shape=0)
# Diagnositc plot, shape zero give good results
par(mfrow=c(2,2))
plot(fit.gum)
estimatorsPOT_0 = fitted(fit.gum)
stdPOT_0 = std.errors(fit.gum)

######################################################################
#return level POT shape NON zero
m10 <- 10*points_year
m100 <- 100*points_year
rl10pot <- th + fit$est[1]/fit$est[2]*((m10*fit$pat)^(fit$est[2])-1)
rl100pot <- th + fit$est[1]/fit$est[2]*((m100*fit$pat)^(fit$est[2])-1)

fit2<-gpd.fit(prod,threshold=th, npy=points_year)
gpd.diag(fit2)

par(mfrow=c(1,2))
#Set good lower and upper bound - manually
gpd.prof(z=fit2, m=10, xlow=2.8e4, xup=3.5e4, npy = points_year, conf = 0.95)
title("Profile Log-likelihood \n of 10-year Return Level")
abline(v=rl10pot)
gpd.prof(z=fit2, m=100, xlow=3.25e4, xup=4.8e4, npy = points_year, conf = 0.95)
title("Profile Log Likelihood \n of 100-year Return Level")
abline(v=rl100pot) #must coincide with MLE

#return level POT shape = 0
rl10_0 <- th + fit.gum$est[1]*log((m10*fit.gum$pat))
rl100_0 <- th + fit.gum$est[1]*log((m100*fit.gum$pat))

###################################################################
# Poisson process
# NOT working with initial threshold th with initial threshold
# information matrix singular

# try higher threshold 0.98 quantile
# keep this threshold, quantile 0.99 not working
th3 <-quantile(prod,0.98)
fit3<-pp.fit(prod,threshold=th3, npy=points_year)
par(mfrow=c(2,2))
pp.diag(fit3)

# compute nb obs above threshold
n <- 0
for (i in c(1:data_size)){
  if (prod[i] > th3){n <- n+1}
}
zeta <- n/data_size

#return level PP shape NON zero
scale_ <- fit3$mle[2] + (th3 - fit3$mle[1])*fit3$mle[3]
rl10pp <- th + scale_/fit3$mle[3]*((m10*zeta)^(fit3$mle[3])-1)
rl100pp <- th + scale_/fit3$mle[3]*((m100*zeta)^(fit3$mle[3])-1)

################################################################
# Bivariate CAPE and SRH
# Monthly maximum
# Max of 8 points per day to 1 point per day and add date sequence 
td <- seq(as.Date("1979/01/01"), as.Date("2015/12/31"), "days")
td <- td[!grepl(x = td, pattern = "-02-29$")]
cape_day <- apply(matrix(cape,ncol=ppd,byrow=TRUE),1,max)
cape_day <- zoo(x = cape_day, order.by = td)
srh_day <- apply(matrix(srh,ncol=ppd,byrow=TRUE),1,max)
srh_day <- zoo(x = srh_day, order.by = td)
# Daily to monthly data taking max per month
cape_m <- daily2monthly(cape_day, FUN=max, na.rm=TRUE) 
srh_m <- daily2monthly(srh_day, FUN=max, na.rm=TRUE) 

# Convert zoo to numeric
cape_day <- as.numeric(cape_day)
cape_m <- as.numeric(cape_m)
srh_day <- as.numeric(srh_day)
srh_m <- as.numeric(srh_m)

# Plot CAPE and SRH
plot(cape_m)
plot(srh_m)

# Fit with symmetric model
fit.mar1 <- fgev(x=cape_m)
fit.mar2 <- fgev(x=srh_m)
res1 <- qgev(pgev(cape_m, loc=fit.mar1$param['loc'],
                  scale=fit.mar1$param['scale'], shape=fit.mar1$param['shape']),1,1,1)
res2 <- qgev(pgev(srh_m, loc=fit.mar2$param['loc'],
                  scale=fit.mar2$param['scale'], shape=fit.mar2$param['shape']),1,1,1)
fbvevd(cbind(res1,res2), cscale=TRUE, cshape=TRUE, cloc=TRUE,
       loc1=1, scale1=1, shape1=1)

cape_srh_m <- cbind(cape_m,srh_m)
fit1 <- fbvevd(cape_srh_m,model="log")
fit1
par(mfrow=c(3,2))
plot(fit1)

#Fit with asymmetric model is singular

# Other bivariate models
fit3 <- fbvevd(cape_srh_m,model="neglog")
par(mfrow=c(3,2))
plot(fit3)
#fit4 <- fbvevd(cape_srh_m,model="bilog") #singular
#fit5 <- fbvevd(cape_srh_m,model="ct")    #singular
#fit6 <- fbvevd(cape_srh_m,model="negbilog") #singular

# AIC
aic1 <- fit1$dev + 2*length(fit1$param)
aic3 <- fit3$dev + 2*length(fit3$param)

# Asymptotic dependence
par(mfrow=c(1,2))
chiplot(cape_srh_m, xlim=c(0.8,1))






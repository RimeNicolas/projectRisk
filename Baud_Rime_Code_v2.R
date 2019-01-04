#install.packages("evd");install.packages("evdbayes");install.packages("coda");install.packages("ismev")
library(evd);library(evdbayes);library(coda);library(ismev)
#install.packages("reliaR");library(reliaR);install.packages("extRemes");library(extRemes)

# Install packages to simplify life
#install.packages("hydroTSM");
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

#Separate by month
prod_monthly <- matrix(prod_m, ncol=12, byrow=TRUE)

# Convert zoo to numeric
x_day <- as.numeric(prod_day)
x_monthly <- as.numeric(prod_monthly)
# Plot the monthly maximum
x_monthly_max <- t(prod_monthly)
matplot(x_monthly_max, lty = 1:12, pch = "o", lend = par("cex"), type="p", col="black")
plot(x_monthly)
# enso monthly
# plot(enso)
enso_monthly_max <- t(matrix(enso, ncol=12, byrow=TRUE))
matplot(enso_monthly_max, lty = 1:12, pch = "o", lend = par("cex"), type="p", col="black")
# Plot Prod monthly accroding to enso_monthly
enso_monthly <- matrix(enso, ncol=12, byrow=TRUE)
matplot( enso_monthly, prod_monthly, pch = "o", lend = par("cex"), type="p", col="black")

for(month in 2:11){
  print(c("SIM NB = ",month))
  x_month = as.numeric(prod_monthly[, month])
  par(mfrow=c(1,1))
  plot(x_month)
  fit_gev = fgev(x_month)
  parameters = fitted(fit_gev) #estimation of the parameters
  print(c("parameters = ",parameters))
  std = std.errors(fit_gev) #std of the parameters (normal dist.)
  print(c("std = ",std))
  par(mfrow=c(2,2))
  plot(fit_gev)
  # Plot profile log likelihood for a better (asymmetric) std approximate
  plot(profile(fit_gev))
}

## For March only, graphs used in the report
years_vec <- c(1979:2015)
x_month = as.numeric(prod_monthly[, 3])

fit_gev = fgev(x_month)
parameters = fitted(fit_gev) #estimation of the parameters
print(c("parameters = ",parameters))
std = std.errors(fit_gev) #std of the parameters (normal dist.)
print(c("std = ",std))

par(mfrow=c(2,2))
plot(fit_gev)

# Plot profile log likelihood for a better (asymmetric) std approximate
plot(profile(fit_gev))
plot(years_vec,x_month, xlab='Years', ylab='Max')
title(main = 'Raw data March')
##

########################################################################
# Check ratio stat with ENSO
for(month in 2:11){
  print(c("SIM NB = ",month))
  x_month = as.numeric(prod_monthly[, month])
  fit_gev = fgev(x_month)
  fit_gev_enso <-fgev(x_month,nsloc=enso[seq(from = month, to = 444, by = 12)])
  ratio <- fit_gev$dev-fit_gev_enso$dev
  p <- 1-pchisq(ratio,1)
  #print(c("ratio = ",ratio))
  print(c("p-value = ",p))
  par(mfrow=c(1,2))
  plot(years_vec,x_month)
  plot(years_vec,enso[seq(from = month, to = 444, by = 12)])
}

# Check ratio stat with TIME

for(month in 2:11){
  print(c("SIM NB = ",month))
  x_month = as.numeric(prod_monthly[, month])
  fit_gev = fgev(x_month)
  fit_gev_time <-fgev(x_month,nsloc=c(1:37))
  ratio <- fit_gev$dev-fit_gev_time$dev
  p <- 1-pchisq(ratio,1)
  #print(c("ratio = ",ratio))
  print(c("p-value = ",p))
  par(mfrow=c(1,2))
  plot(years_vec,x_month)
  plot(years_vec,enso[seq(from = month, to = 444, by = 12)])
}

# Check correlation with ENSO
for(month in 1:12){
  print(c("SIM NB = ",month))
  x_month <- as.numeric(prod_monthly[, month])
  mat1 <- cbind(x_month,enso[seq(from = month, to = 444, by = 12)])
  plot(mat1)
  corr <- cor(mat1)
  print(c("Correlation = ",corr[2]))
  #mat2 <- cbind(x_month,abs(enso[seq(from = month, to = 444, by = 12)]))
  #plot(mat2)
  #corr_abs <- cor(mat2)
  #print(c("Correlation ABSOLU= ",corr_abs))
}

#######################################################
#MCMC
# choose month between 1 and 12
month <- 10
x_m = as.numeric(prod_monthly[, month])
init <- c(5000,1000,0)
mat = diag(c(1e8,1e5,1e5))
pn = prior.norm(mean=c(0,0,0),cov=mat) # normal distribution, flat prior (uninformative)

if (month == 1){psd = c(400,0.3,0.5)}# sd for proposal dist (normal dist)
if (month == 2){psd = c(1000,0.4,0.5)}
if (month == 3){psd = c(1800,0.4,0.5)}
if (month == 4){psd = c(2900,0.4,0.3)}
if (month == 5){psd = c(2900,0.4,0.3)}
if (month == 6){psd = c(2200,0.3,0.3)}
if (month == 7){psd = c(2000,0.4,0.4)}
if (month == 8){psd = c(1500,0.4,0.3)}
if (month == 9){psd = c(1800,0.4,0.3)}
if (month == 10){psd = c(1200,0.4,0.5)}
if (month == 11){psd = c(1200,0.4,0.4)}
if (month == 12){psd = c(400,0.4,0.5)}
 
# initial value, prior distribution: normal, likelihood: gev, psd: 
post = posterior(5000,init=init,prior=pn,lh="gev",data=x_m,psd=psd)

mc = mcmc(post)
plot(mc)
attr(mc,"ar")

MCMC2<-mcmc(post[-c(1:500),])
acf(MCMC2)

MCMC3<-mcmc(post[500 + 9*c(1:500),])
acf(MCMC3)
apply(MCMC3,2,mean)
apply(MCMC3,2,sd)

# Return level MCMC
hist(as.numeric(MCMC3[,3]),nclass=20,prob=T,main="Histogram of xi",xlab="xi")
u.50<-mc.quant(MCMC3,p=0.98,lh="gev")
u.100<-mc.quant(MCMC3,p=0.99,lh="gev")
xlim_ <- range(c(0:5e4))
par(mfrow=c(1,2))
hist(u.50,nclass=200,prob=T,xlab="50-year return level",xlim=xlim_)
hist(u.100,nclass=400,prob=T,xlab="100-year return level",xlim=xlim_)

###############################################
# r largest statistics per month
r_ <- 1
jan <- matrix(nrow = years, ncol = r_)
feb <- matrix(nrow = years, ncol = r_)
mar <- matrix(nrow = years, ncol = r_)
apr <- matrix(nrow = years, ncol = r_)
may <- matrix(nrow = years, ncol = r_)
jun <- matrix(nrow = years, ncol = r_)
jul <- matrix(nrow = years, ncol = r_)
aug <- matrix(nrow = years, ncol = r_)
sep <- matrix(nrow = years, ncol = r_)
oct <- matrix(nrow = years, ncol = r_)
nov <- matrix(nrow = years, ncol = r_)
dec <- matrix(nrow = years, ncol = r_)

index_ <- 1
for (i in c(1:37)){
  jan[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd 
  
  feb[i,] <- rev(tail(sort(prod[index_:(index_+28*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 28*ppd
  
  mar[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd
  
  apr[i,] <- rev(tail(sort(prod[index_:(index_+30*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 30*ppd
  
  may[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd
  
  jun[i,] <- rev(tail(sort(prod[index_:(index_+30*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 30*ppd
  
  jul[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd
  
  aug[i,] <- rev(tail(sort(prod[index_:(index_+30*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd
  
  sep[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 30*ppd
  
  oct[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd
  
  nov[i,] <- rev(tail(sort(prod[index_:(index_+30*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 30*ppd
  
  dec[i,] <- rev(tail(sort(prod[index_:(index_+31*ppd)], decreasing=FALSE),r_))
  index_ <- index_ + 31*ppd

}
print(index_)

# Plot months
mon <- feb
par(mfrow=c(1,1))
time_ <- rep(1:years,each=r_)
plot(time_,t(mon))

# Fit r-largest stat 
fit0<-rlarg.fit(mon) # constant
fit1<-rlarg.fit(mon,ydat=matrix(time_,ncol=r_),mul=c(1)) # linear

fit0$mle
fit1$mle


########################################################################
# Peaks over threshold
jan <- c(1:(ppd*31*37))
feb <- c(1:(ppd*28*37))
mar <- c(1:(ppd*31*37))
apr <- c(1:(ppd*30*37))
may <- c(1:(ppd*31*37))
jun <- c(1:(ppd*30*37))
jul <- c(1:(ppd*31*37))
aug <- c(1:(ppd*31*37))
sep <- c(1:(ppd*30*37))
oct <- c(1:(ppd*31*37))
nov <- c(1:(ppd*30*37))
dec <- c(1:(ppd*31*37))

index_ <- 1
for (i in c(1:37)){
  jan[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
  feb[(1+(28*ppd*(i-1))):(28*ppd*i)] <- prod[index_:(index_+28*ppd-1)]
  index_ <- index_ + 28*ppd 
  
  mar[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
  apr[(1+(30*ppd*(i-1))):(30*ppd*i)] <- prod[index_:(index_+30*ppd-1)]
  index_ <- index_ + 30*ppd 
  
  may[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
  jun[(1+(30*ppd*(i-1))):(30*ppd*i)] <- prod[index_:(index_+30*ppd-1)]
  index_ <- index_ + 30*ppd 
  
  jul[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
  aug[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
  sep[(1+(30*ppd*(i-1))):(30*ppd*i)] <- prod[index_:(index_+30*ppd-1)]
  index_ <- index_ + 30*ppd 
  
  oct[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
  nov[(1+(30*ppd*(i-1))):(30*ppd*i)] <- prod[index_:(index_+30*ppd-1)]
  index_ <- index_ + 30*ppd 
  
  dec[(1+(31*ppd*(i-1))):(31*ppd*i)] <- prod[index_:(index_+31*ppd-1)]
  index_ <- index_ + 31*ppd 
  
}
print(index_)

# Choose month
mon <- dec; month <- 12

par(mfrow=c(1,1))
plot(mon)
qu.min <- quantile(mon, 0.5) # median value
qu.max <- quantile(mon,(length(mon)-30)/length(mon))
print(paste0("median: ",qu.min,", quantile max: ",qu.max))
mrlplot(mon, tlim=c(qu.min, qu.max))
par(mfrow=c(1,2))
tcplot(mon,tlim=c(qu.min, qu.max))

# Choose threshold by hand giving a good variance bias trade-off
if (month == 1){th <- 1.2e3; points_month <- 8*31; extr_interval <- c(500,3000)}
if (month == 2){th <- 4e3; points_month <- 8*28; extr_interval <- c(200,4700)}
if (month == 3){th <- 9e3; points_month <- 8*31; extr_interval <- c(4e3,11e3)}
if (month == 4){th <- 12e3; points_month <- 8*30; extr_interval <- c(4e3,12e3)}
if (month == 5){th <- quantile(mon,0.9); points_month <- 8*31; extr_interval <- c(3e3,12e3)}
if (month == 6){th <- 10e3; points_month <- 8*30; extr_interval <- c(4e3,12e3)}
if (month == 7){th <- 12e3; points_month <- 8*31; extr_interval <- c(9e3,14e3)}
if (month == 8){th <- 8e3; points_month <- 8*31; extr_interval <- c(6e3,12e3)}
if (month == 9){th <- 10e3; points_month <- 8*30; extr_interval <- c(8e3,14e3)}
if (month == 10){th <- 8e3; points_month <- 8*31; extr_interval <- c(6e3,11e3)}
if (month == 11){th <- 5e3; points_month <- 8*30; extr_interval <- c(3e3,6e3)}
if (month == 12){th <- 2e3; points_month <- 8*31; extr_interval <- c(1e3,2.4e3)}

fit<-fpot(mon,model="gpd",threshold=th,npp=points_month)
par(mfrow=c(2,2))
plot(fit)

estimatorsPOT = fitted(fit) #estimation of the parameters
stdPOT = std.errors(fit) #std of the parameters (normal dist.)

# Profile log likelihood 
par(mfrow=c(1,2))
plot(profile(fit))
abline(v=0,col=2,lty=2)
# Zero shape should be considered
#fit.gum<-fpot(mon, threshold=th, npp=points_month, shape=0)
# Diagnositc plot, shape zero
#par(mfrow=c(2,2))
#plot(fit.gum)
#estimatorsPOT_0 = fitted(fit.gum)
#stdPOT_0 = std.errors(fit.gum)

print(c("Threshold = ",th))
print(c("estimators POT = ",estimatorsPOT))
print(c("std POT = ",stdPOT))
#print(c("estimators POT shape:0 = ",estimatorsPOT_0))
#print(c("std POT shape: 0 = ",stdPOT_0))

######################################################################
# Extremal index 
par(mfrow=c(1,1))
exiplot(jan, tlim=extr_interval)
exiplot(jan, tlim=extr_interval, r=0, add=T, lty=2)

######################################################################
#return level POT shape NON zero
m50 <- 50*points_month
m100 <- 100*points_month
rl50pot <- th + fit$est[1]/fit$est[2]*((m50*fit$pat)^(fit$est[2])-1)
rl100pot <- th + fit$est[1]/fit$est[2]*((m100*fit$pat)^(fit$est[2])-1)

#return level POT shape = 0
#rl50_0 <- th + fit.gum$est[1]*log((m50*fit.gum$pat))
#rl100_0 <- th + fit.gum$est[1]*log((m100*fit.gum$pat))

###################################################################
# Poisson process
th3 <- th
#th3 <-quantile(mon,0.97)
fit3<-pp.fit(mon,threshold=th3, npy=points_month)
par(mfrow=c(2,2))
pp.diag(fit3)

# compute nb obs above threshold
n <- 0
for (i in c(1:length(mon))){
  if (mon[i] > th3){n <- n+1}
}
zeta <- n/length(mon)

#return level PP shape NON zero
scale_ <- fit3$mle[2] + (th3 - fit3$mle[1])*fit3$mle[3]
rl50pp <- th3 + scale_/fit3$mle[3]*((m50*zeta)^(fit3$mle[3])-1)
rl100pp <- th3 + scale_/fit3$mle[3]*((m100*zeta)^(fit3$mle[3])-1)

print(c("Threshold POT = ",th))
print(c("Threshold PP = ",th3))
print(c("rl50 POT = ",rl50pot))
print(c("rl100 POT = ",rl100pot))
#print(c("rl50 POT shape:0 = ",rl50_0))
#print(c("rl100 POT shape:0 = ",rl100_0))
print(c("rl50 PP = ",rl50pp))
print(c("rl100 PP = ",rl100pp))

# Fig for report
par(mfrow=c(2,2))
mrlplot(mon, tlim=c(qu.min, qu.max),main=NULL)
exiplot(jan, tlim=extr_interval)
exiplot(jan, tlim=extr_interval, r=0, add=T, lty=2)
tcplot(mon,tlim=c(qu.min, qu.max))

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
par(mfrow=c(1,2))
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

########################################################################################
# Bivariate per month
cape_monthly <- matrix(cape_m, ncol=12, byrow=TRUE)
srh_monthly <- matrix(srh_m, ncol=12, byrow=TRUE)

# month number
month <- 2
cape_month <- as.numeric(cape_monthly[, month])
srh_month <- as.numeric(srh_monthly[, month])

# Plot CAPE and SRH
par(mfrow=c(1,2))
plot(cape_month)
plot(srh_month)

# Fit with symmetric model
fit.mar1 <- fgev(x=cape_month)
fit.mar2 <- fgev(x=srh_month)
res1 <- qgev(pgev(cape_month, loc=fit.mar1$param['loc'],
                  scale=fit.mar1$param['scale'], shape=fit.mar1$param['shape']),1,1,1)
res2 <- qgev(pgev(srh_month, loc=fit.mar2$param['loc'],
                  scale=fit.mar2$param['scale'], shape=fit.mar2$param['shape']),1,1,1)
fbvevd(cbind(res1,res2), cscale=TRUE, cshape=TRUE, cloc=TRUE,
       loc1=1, scale1=1, shape1=1)

cape_srh_month <- cbind(cape_month,srh_month)
fit1 <- fbvevd(cape_srh_month,model="log")
fit1
par(mfrow=c(3,2))
plot(fit1)

fit2 <- fbvevd(cape_srh_month,model="alog")
fit2
par(mfrow=c(3,2))
plot(fit2)

# Other bivariate models
fit3 <- fbvevd(cape_srh_month,model="neglog")
par(mfrow=c(3,2))
plot(fit3)
#fit4 <- fbvevd(cape_srh_month,model="bilog") #singular
fit5 <- fbvevd(cape_srh_month,model="ct")    #singular
#fit6 <- fbvevd(cape_srh_month,model="negbilog") #singular

# AIC
aic1 <- fit1$dev + 2*length(fit1$param)
aic2 <- fit2$dev + 2*length(fit2$param)
aic3 <- fit3$dev + 2*length(fit3$param)
aic5 <- fit5$dev + 2*length(fit5$param)

# Asymptotic dependence
par(mfrow=c(1,2))
chiplot(cape_srh_month, xlim=c(0.5,1))

########################################################################################
# log PROD
for(month in 2:11){
  print(c("SIM NB = ",month))
  x_month = as.numeric(log(prod_monthly[, month]))
  par(mfrow=c(1,1))
  plot(x_month)
  fit_gev = fgev(x_month)
  parameters = fitted(fit_gev) #estimation of the parameters
  print(c("parameters = ",parameters))
  std = std.errors(fit_gev) #std of the parameters (normal dist.)
  print(c("std = ",std))
  par(mfrow=c(2,2))
  plot(fit_gev)
  # Plot profile log likelihood for a better (asymmetric) std approximate
  plot(profile(fit_gev))
}




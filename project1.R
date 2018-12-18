install.packages("evd");install.packages("evdbayes");install.packages("coda")
library(evd);library(evdbayes);library(coda)
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
days_y <- 365
years <- 37
data_size <- ppd*days*years

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
x_m <- as.numeric(prod_m)

# Plot the monthly maximum
plot(x_m)

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

#MCMC
init <- c(5000,1000,0)
mat = diag(c(1e11,1e8,1e8))
pn = prior.norm(mean=c(0,0,0),cov=mat)

psd = c(800,0.10,0.12) #variance of estimators
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



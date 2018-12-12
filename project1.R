#install.packages("evd");install.packages("evdbayes");
library(evd);library(evdbayes);
#install.packages("reliaR");library(reliaR)
#install.packages("extRemes");library(extRemes)

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

# Convert zoo to numeric
x_m <- as.numeric(prod_m)

# QQplot
months <- length(x_m)
qqplot(qgumbel(c(1:months)/(months+1)),x_m)
qqline(x_m,distribution=qgumbel)

# Fit GEV 
fit_gev = fgev(x_m)
estimators = fitted(fit_gev)
std = std.errors(fit_gev)
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(fit_gev)
plot(profile(fit_gev))

plot(profile(fgev(m1,prob=1-exp(-0.1)),"quantile"))
qgumbel(exp(-0.1),loc=2.82,scale=0.31)




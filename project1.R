#install.packages("evd");install.packages("evdbayes");
library(evd);library(evdbayes);
#install.packages("reliaR");library(reliaR)
#install.packages("extRemes");library(extRemes)

years <- 50
days <- 365
data_size <- years*days
x1 <- rnorm(data_size)
x2 <- 0.99*rnorm(data_size,5,1) + 0.01*rnorm(data_size,5,100)
x3 <- rcauchy(data_size)
x4 <- rbeta(data_size,6,9)
x5 <- 1 - 1/log(runif(data_size))
x6 <- log(x5)
x  <- cbind(x1,x2,x3,x4,x5,x6)
nb_dist = length(x[1,])

m <- matrix(nrow=nb_dist,ncol=years)

for (i in 1:nb_dist){
  m[i,] <- apply(matrix(x[,i],ncol=days),1,max)
}
#From here consider only normal distribution dist_test
dist_test = 1
m1 = m[dist_test,]
qqplot(qgumbel(c(1:years)/(years+1)),m[dist_test,])
qqline(m[dist_test,],distribution=qgumbel)

f_gev = fgev(m[dist_test,])
estimators = fitted(f_gev)
mystd = std.errors(f_gev)
#par(mfrow=c(2,2))
par(mfrow=c(1,1))
plot(f_gev)
plot(profile(f_gev))


fit1 = fgev(m1)
fit2 = fgev(m1,shape=0)
ratio<-fit2$dev-fit1$dev
qchisq(0.95,1)
fit1
fit2

f_gev
plot(profile(fgev(m1,prob=1-exp(-0.1)),"quantile"))
qgumbel(exp(-0.1),loc=2.82,scale=0.31)

#d
x7<-matrix(rweibull(10^7,2,1), nrow=10000)
x8<-matrix(rweibull(10^7,2,1), nrow=10000)^2

m7 <- apply(x7, 1, max)
m8 <- apply(x8, 1, max)

par(mfrow=c(1,1))
N<-seq(from=10,to=10000,by=100)
xi.1<-c()
xi.2<-c()
for(i in 1:length(N)){
  fit.1<-fgev(m7[1:N[i]])
  xi.1[i]<-fit.1$est[3]
  fit.2<-fgev(m8[1:N[i]])
  xi.2[i]<-fit.2$est[3]
}
plot(N,xi.1,type="l",ylim=c(min(xi.1,xi.2,-0.5/log(10)),max(xi.1,xi.2)),
     ylab="shape parameter")
points(N,xi.2,type="l",col="red")
abline(h=0,col="blue")
points(N,-0.5/log(N),type="l",col="brown")


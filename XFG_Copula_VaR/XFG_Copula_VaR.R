## 1
# clear all variables
rm(list = ls(all = TRUE))
graphics.off()

#Set working directory
setwd("~/Desktop/IRTG/Teaching/R")

# install and load packages
libraries = c("copula","PerformanceAnalytics")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

## 2
#Read data
names = c("WFC","JPM","BAC","C","BK","STT","GS","MS")
x0 = as.matrix(read.csv(file = "Returns.csv"))
date = read.csv(file = "Date.csv",colClasses = "Date")

## 3
#Get pseudo observations (~U[0,1])
u = pobs(x0)
gauss_fit = fitCopula(normalCopula(dim=ncol(x0)),u,method="mpl")
t_fit     = fitCopula(tCopula(dim=ncol(x0)),u,method="mpl")
clay_fit  = fitCopula(claytonCopula(dim=ncol(x0)),u,method="mpl")
gumb_fit  = fitCopula(gumbelCopula(dim=ncol(x0)),u,method="mpl")

## 4
#Simulate values
N = 1000

gauss_sim = rCopula(N,normalCopula(coef(gauss_fit),dim=ncol(x0)))
t_sim     = rCopula(N,tCopula(coef(t_fit),dim=ncol(x0)))
clay_sim  = rCopula(N,claytonCopula(coef(clay_fit),dim=ncol(x0)))
gumb_sim  = rCopula(N,gumbelCopula(coef(gumb_fit),dim=ncol(x0)))

## 5
#Correlation matrices
chart.Correlation(gauss_sim)
chart.Correlation(t_sim)
chart.Correlation(clay_sim)
chart.Correlation(gumb_sim)

pairs(gauss_sim,pch=19,cex=0.01)
pairs(t_sim,pch=19,cex=0.01)
pairs(clay_sim,pch=19,cex=0.01)
pairs(gumb_sim,pch=19,cex=0.01)

## 6
#Simulated pobs + empirical quantiles = meta distribution
gauss_dis = t_dis = clay_dis = gumb_dis = matrix(0, nrow=N, ncol=ncol(x0))
for (j in 1:ncol(x0)){
  gauss_dis[,j] = quantile(x0[,j],gauss_sim[,j])
  #t_dis[,j]     = quantile(x0[,j],t_sim[,j])
  clay_dis[,j]  = quantile(x0[,j],clay_sim[,j])
  gumb_dis[,j]  = quantile(x0[,j],gumb_sim[,j])
}

## 7
# Construct distribution of equally weighted stock portfolio

gauss_pf1 = apply(gauss_dis,1,mean)
t_pf1     = apply(t_dis,1,mean)
clay_pf1  = apply(clay_dis,1,mean)
gumb_pf1  = apply(gumb_dis,1,mean)

## 8
# Density plots + VaR and ES
tau = 0.05

plot(density(gauss_pf1,bw=0.01),main="Gaussian Copula",xlim=c(-0.15,0.1))
abline(v=quantile(gauss_pf1,tau),col="blue")
abline(v=mean(gauss_pf1[gauss_pf1<quantile(gauss_pf1,tau)]),col="red")

#plot(density(t_pf1,bw=0.01),main="t Copula",xlim=c(-0.15,0.1))
#abline(v=quantile(t_pf1,tau),col="blue")
#abline(v=mean(t_pf1[t_pf1<quantile(t_pf1,tau)]),col="red")

plot(density(clay_pf1,bw=0.01),main="Clayton Copula",xlim=c(-0.15,0.1))
abline(v=quantile(clay_pf1,tau),col="blue")
abline(v=mean(clay_pf1[clay_pf1<quantile(clay_pf1,tau)]),col="red")

plot(density(gumb_pf1,bw=0.01),main="Gumbel Copula",xlim=c(-0.15,0.1))
abline(v=quantile(gumb_pf1,tau),col="blue")
abline(v=mean(gumb_pf1[gumb_pf1<quantile(gumb_pf1,tau)]),col="red")

#Historical
plot(density(apply(x0,1,mean),bw=0.01),main="Historical Simulation",xlim=c(-0.15,0.1))
abline(v=quantile(apply(x0,1,mean),tau),col="blue")
abline(v=mean(apply(x0,1,mean)[apply(x0,1,mean)<quantile(apply(x0,1,mean),tau)]),col="red")


## 9
#GARCH(1,1)

## 10
#Dynamic Setting

## 11
#Portfolio with nonlinear exposure
delta = c(5,5,-5,-5,10,10,-10,-10)
gamma = rep(-0.8,8)

gauss_pf2 = gauss_dis %*% delta + 1/2 * gauss_dis^2 %*% gamma
t_pf2     = t_dis %*% delta + 1/2 * t_dis^2 %*% gamma
clay_pf2  = clay_dis %*% delta + 1/2 * clay_dis^2 %*% gamma
gumb_pf2  = gumb_dis %*% delta + 1/2 * gumb_dis^2 %*% gamma

#Gauss
plot(density(gauss_pf2,bw=0.05),main="Gaussian Copula - Pf2",xlim=c(-2,2))
abline(v=quantile(gauss_pf2,tau),col="blue")
abline(v=mean(gauss_pf2[gauss_pf2<quantile(gauss_pf2,tau)]),col="red")

#Clayton
plot(density(clay_pf2,bw=0.05),main="Clayton Copula - Pf2",xlim=c(-2,2))
abline(v=quantile(clay_pf2,tau),col="blue")
abline(v=mean(clay_pf2[clay_pf2<quantile(clay_pf2,tau)]),col="red")

#Gumbel
plot(density(gumb_pf2,bw=0.05),main="Gumbel Copula - Pf2",xlim=c(-2,2))
abline(v=quantile(gumb_pf2,tau),col="blue")
abline(v=mean(gumb_pf2[gumb_pf2<quantile(gumb_pf2,tau)]),col="red")

#Historical
hist_pf2  = x0 %*% delta + 1/2 * x0^2 %*% gamma

plot(density(hist_pf2,bw=0.05),main="Historical Simulation - Pf2",xlim=c(-2,2))
abline(v=quantile(hist_pf2,tau),col="blue")
abline(v=mean(hist_pf2[hist_pf2<quantile(hist_pf2,tau)]),col="red")

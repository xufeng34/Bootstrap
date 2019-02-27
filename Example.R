rm(list=ls())
# install.packages(c('forecast','TSA','qcc','nortest'))
library(forecast)
library(TSA)
library(qcc)
library(nortest)
# setwd("D:/")
source("Function.R")
options(scipen = 200) #不采用科学计数法
dat = read.table('绝缘体.txt') 
X = dat$V1
X 
n=length(X);n
par(mar=c(4,4,1,1)+0.5)
plot(X, type='o',pch=16, ylab = '')
# hist(X, freq = F, yaxt="n", main='',xlab = '',ylab = '' )
hist(X, freq = FALSE, breaks = 20, main='',xlab = '',ylab = '' )
# axis(side = 2,at = seq(0,10,by=1),labels = seq(0,1,by=0.1) )
lines(density(X))
abline(v = c(2800,(2800+5800)/2,5800),lty = 1);abline(v = mean(X),lty = 2);
text(2900,0.0009,'LCL',xpd = T,cex=0.75);text(5700,0.0009,'UCL',xpd = T,cex=0.75);
text((2800+5800)/2-50,0.0009,'M',xpd = T,cex=0.75);text(mean(X)-50,0.0009,expression(mu),xpd = T,cex=0.75);
box()
# 3 检验是否存在自相关----
## 自相关
par(mar=c(4,4,1,1)+0.5)
Acf(X)  #600*400/6.24*4.18PDF
## 偏自相关
pacf(X)


## 利用残差图判断
# fitmodel <- auto.arima(X);fitmodel
fitmodel <- arima(X,order=c(1,0,0));fitmodel
# names(fitmodel)
residuals0 <- fitmodel$residuals

obj <- qcc(residuals0, type="xbar.one",plot = TRUE)
CL <- obj$center
LCL <- obj$limits[1]
UCL <- obj$limits[2]

par(mar=c(4,2,1,2)+0.5)
plot(residuals0,type='o',pch = 16,xlab = "")
points(which(residuals0>UCL|residuals0<LCL),residuals0[residuals0>UCL|residuals0<LCL],pch = 16, col = 'red')
abline(h = CL);abline(h = c(LCL,UCL),lty = 2)
text(220,LCL,'LCL',xpd = T,cex=0.75);text(220,UCL,'UCL',xpd = T,cex=0.75);text(220,CL,'CL',xpd = T,cex=0.75)

## 去掉超出控制限的点后重新构造控制图
flag = 0
while(length(flag)!=0){
  X <- X[-which(residuals0>UCL|residuals0<LCL)]
  fitmodel <- arima(X,order=c(1,0,0));fitmodel
  residuals0 <- fitmodel$residuals
  
  obj <- qcc(residuals0, type="xbar.one",plot = TRUE)
  CL <- obj$center
  LCL <- obj$limits[1]
  UCL <- obj$limits[2]
  
  flag = which(residuals0>UCL|residuals0<LCL)
}
par(mar=c(4,2,1,2)+0.5)
plot(residuals0,type='o',pch = 16,ylim = c(-1100,1100),xlab = "")
points(which(residuals0>UCL|residuals0<LCL),residuals0[residuals0>UCL|residuals0<LCL],pch = 16, col = 'red')
abline(h = CL);abline(h = c(LCL,UCL),lty = 2)
text(215,LCL,'LCL',xpd = T,cex=0.75);text(215,UCL,'UCL',xpd = T,cex=0.75);text(215,CL,'CL',xpd = T,cex=0.75)

fitmodel <- arima(X,order=c(1,0,0));fitmodel

### interval  ----
U = 5800; L = 2800; 

fitm <- arima(X, order = c(1,0,0), include.mean = TRUE)
PCIB <- BootPCI(X, fitm = fitm, B = 1000, U, L, n = n)

# PCIhat <- Cp(mu = mean(X), sigma = sd(X), U=U, L=L)
phihat <- fitm$coef[1]
sigmahat = trsigma(x = X)/sqrt(1-phihat)
PCIhat <- Cpk(mu = fitm$coef[2], sigma = sigmahat, U = U, L = L)

  
Result <- BootstrapInterval(BootSamp = PCIB$PCIB, PCIHat = PCIhat, alpha = 0.05)
Result 





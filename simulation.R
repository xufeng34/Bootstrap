rm(list=ls())
# install.packages(c("TSA","doParallel"))
library(TSA)
library(doParallel)
# library(foreach)
# setwd("D:/")
source("Function.R")

time <- Sys.time()
# U = 61; L = 40; mu = 50; sigma = 2; phi = 0.4; n = 100;
# distribution = "chisq"; para1 = 20; para2 = NULL; tau = 100
# distribution = "F"; para1 = 30; para2 = 12; tau = 100
# distribution = "norm"; para1 = 0; para2 = 1; tau = 100
# distribution = "gamma"; para1 = 4; para2 = 2; tau = 100

simulation <- function(U, L, mu, sigma, phi, n, distribution, para1, para2, tau = 100){
  SB_perc = 0;  
  SB_width = vector();  
  PCI.Value = Cpk(mu, sigma, U, L)
  iteration = 1000
  kk = 0
  while(kk < iteration){
    kk = kk + 1
    if(distribution == "F") dat <- GenerateData(mu, sigma, phi, distribution = "F", para1, para2, n, tau)
    if(distribution == "gamma") dat <- GenerateData(mu, sigma, phi, distribution = "gamma", para1, para2, n, tau)
    if(distribution == "chisq") dat <- GenerateData(mu, sigma, phi, distribution = "chisq", para1, para2, n, tau)
    if(distribution == "norm") dat <- GenerateData(mu, sigma, phi, distribution = "norm", para1, para2, n, tau)
    
    fitm <- tryCatch(arima(dat, order = c(1,0,0), include.mean = TRUE), error = function(e) e)
    if(length(class(fitm)) > 1){ if(class(fitm)[2] == "error") {kk = kk-1; next} }
    
    phihat <- fitm$coef[1]
    sigmahat = trsigma(x = dat)/sqrt(1-phihat)
    
    PCIhat <- Cpk(mu = fitm$coef[2], sigma = sigmahat, U = U, L = L)
    
    PCIB <- BootPCI(dat, fitm = fitm, B = 1000, U, L, n = n)
    
    Result <- BootstrapInterval(BootSamp = PCIB$PCIB, PCIHat = PCIhat, alpha = 0.05)
    
    if(PCI.Value >= Result[1] & PCI.Value <= Result[2]) SB_perc <- SB_perc + 1
    SB_width[kk] <- Result[3]; 
  }
  
  SB_width = SB_width[SB_width != 0]; nSB = length(SB_width)
  int <- data.frame("Coverage" = SB_perc/kk,"width" = sum(SB_width)/nSB, "SD" = sd(SB_width))
  int <- round(int,3)
  int <- data.frame("distr"=distribution, 'n' = n,"phi" = phi,"mu"=mu,"sigma"=sigma,int)
  int
}

detectCores()         
cl <- makeCluster(5)  
registerDoParallel(cl)


time <- Sys.time()
case1 <- foreach(phi = rep(c(-0.8,-0.4,0,0.4,0.8),4),sigma = c(rep(2,10),rep(3.5,10)),mu = c(rep(50,5),rep(52,5),rep(50,5),rep(52,5)),  .combine='rbind',.packages=c("TSA",'doParallel')) %dopar%
  simulation(U = 61, L = 40, mu = mu, sigma = sigma, phi = phi, distribution = "norm", para1 = 0, para2 = 1, n = 100, tau = 100)
write.table(case1, paste0("D:/norm_n100.csv"), append = TRUE,sep = ',',row.names = FALSE)
f(time)
Sys.time()

time <- Sys.time()
case2 <- foreach(phi = rep(c(-0.8,-0.4,0,0.4,0.8),4),sigma = c(rep(2,10),rep(3.5,10)),mu = c(rep(50,5),rep(52,5),rep(50,5),rep(52,5)), .combine='rbind',.packages=c("TSA",'doParallel')) %dopar%
  simulation(U = 61, L = 40, mu = mu, sigma = sigma, phi = phi, distribution = "chisq", para1 = 20, para2 = NULL, n = 100, tau = 100)
write.table(case2, paste0("D:/chisq_n100.csv"), append = TRUE,sep = ',',row.names = FALSE)
f(time)
Sys.time()

time <- Sys.time()
case3 <- foreach(phi = rep(c(-0.8,-0.4,0,0.4,0.8),4),sigma = c(rep(2,10),rep(3.5,10)),mu = c(rep(50,5),rep(52,5),rep(50,5),rep(52,5)), .combine='rbind',.packages=c("TSA",'doParallel')) %dopar%
  simulation(U = 61, L = 40, mu = mu, sigma = sigma, phi = phi, distribution = "gamma", para1 = 4, para2 = 2, n = 100, tau = 100)
write.table(case3, paste0("D:/gamma_n100.csv"), append = TRUE,sep = ',',row.names = FALSE)
f(time)
Sys.time()

stopCluster(cl)
f(time)
Sys.time()
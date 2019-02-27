library(TSA)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
#
BootstrapInterval <- function(BootSamp, PCIHat, alpha = 0.05){
  # SB
  Za = qnorm(1-alpha/2)
  S <- sd(BootSamp)
  L_SB <- PCIHat - Za*S
  U_SB <- PCIHat + Za*S
  width_SB <- U_SB-L_SB
  SB <- c(L_SB, U_SB, width_SB)
  names(SB) <- c("Lower","Upper","width")
  SB
  
}
#
Cpk <- function(mu, sigma, U, L){
  min(U-mu,mu-L)/(3*sigma)
}
# generate data ----
# U = 61; L = 40; 
# mu = 50; sigma = 2; phi = 0.4; distribution = "F"; para1 = 30; para2 = 12; n = 100; tau = 100
GenerateData <- function(mu, sigma, phi, distribution, para1, para2, n, tau){
  # Gamma ----
  # rho = phi = 0.4
  # tau = 100; n = 100
  # para1 = 4; para2 = 2
  # mu = 50; sigma = 2
  # M = para1/para2; S = para1/(para2^2)
  rho = phi
  if(distribution == "F"){
    M = para2/(para2-2);
    S = sqrt((2*para2^2*(para1+para2-2))/(para1*(para2-2)^2*(para2-4)));
  }
  if(distribution == "gamma"){
    M = para1/para2
    S = sqrt(para1/(para2^2))
  }
  if(is.null(para2)|distribution == "chisq"){
    M = para1;
    S = sqrt(2*para1);  
  }
  if(distribution == "norm"){
    M = para1
    S = para2
  }

  flag = TRUE
  while(flag){
    maxtimes = 0
    while(flag&maxtimes<2000){
      m = 100
      rhohat = vector()
      for(i in 1:m){
        w = arima.sim(model = list(ar = phi), n = n + tau, sd = sqrt(1-phi^2))
        z = w[tau+(1:n)]
        u = pnorm(z)
        if(distribution == "F"){y = qf(u,para1,para2) }
        if(distribution == "gamma"){y = qgamma(u,para1,para2)  }
        if(is.null(para2)|distribution == "chisq"){    y = qchisq(u,para1)    }
        if(distribution == "norm"){    y = qnorm(u,para1,para2)  }
        x = ((y-M)/S)*sigma+mu
        rhohat[i] <- acf(x, plot = FALSE)[[1]][1] # bear in mind using package TSA
      }
      rhobar <- mean(rhohat)
      serhobar <- sqrt(sum((rhohat-rhobar)^2)/(m*(m+1)))
      Delta <- abs(rhobar-rho)
      flag = !(Delta <= serhobar);flag
      if(flag) phi = ifelse(rhobar>rho,phi-Delta,phi+Delta); 
      maxtimes = maxtimes + 1
      # print(c(Delta,serhobar))
      # print(paste0("k0=",maxtimes))
      # print(phi)
    }
  }
  # w = arima.sim(model = list(ar = phi), n = n + tau, sd = sqrt(1-phi^2))
  # z = w[tau+(1:n)]
  # u = pnorm(z)
  # if(distribution == "F"){y = qf(u,para1,para2) }
  # if(distribution == "gamma"){y = qgamma(u,para1,para2)  }
  # if(is.null(para2)|distribution == "chisq"){    y = qchisq(u,para1)    }
  # if(distribution == "norm"){ y = qnorm(u,para1,para2)  }
  # x = ((y-M)/S)*sigma+mu
  # x
  return(x)
}
# GenerateData(mu, sigma, phi, distribution, para1, para2, n, tau)
BootPCI <- function(dat, fitm = NULL, B = 1000,U,L, n ){
  if(is.null(fitm)){
    n = length(dat)
    fitm <- tryCatch(arima(dat, order = c(1,0,0), include.mean = TRUE),error = function(e) e)
    if(length(class(fitm)) > 1){
      if(class(fitm)[2] == "error") return(cat(paste0(x))) 
    }
  }
  
  phihat <- fitm$coef[1]; 
  intercept <- fitm$coef[2]
  e <- fitm$residuals
  st_e <- scale(e, center = TRUE, scale = FALSE) # standard residuals
  ce <- mean(e)
  PCIB <- vector()
  k = 0
  while(k < B){
    k = k + 1
    re_st_e <- sample(st_e, size = n, replace = TRUE)
    non_ce <- re_st_e + ce
    X0 <- sample(dat,1)
    X <- vector()
    for(i in 1:n){
      if(i == 1){
        X[i] <- phihat*X0 + non_ce[i] + (1-phihat)*intercept
      }else{
        X[i] <- phihat*X[i-1] + non_ce[i] + (1-phihat)*intercept
      }
    }
    muB = mean(X)
    
    fitm2 <- tryCatch(arima(X, order = c(1,0,0), include.mean = TRUE),error = function(e) e)
    if(length(class(fitm2)) > 1){ if(class(fitm2)[2] == "error") {k = k - 1; next } }
    phihat2 <- fitm2$coef[1]
    sigmaB = trsigma(x = X)/sqrt(1-phihat2)
    # sigmaB = sd(X)
    PCIB[k] <- Cpk(muB,sigmaB,U,L)
  }
  PCIBbar <- mean(PCIB)
  return(list("PCIB" = PCIB, "PCIBbar" = PCIBbar))
}
#

trsigma <- function(x){
  n = length(x)
  Rbar <- sum(abs(diff(x)))/(n-1)
  sigmahat <- Rbar/(2/sqrt(pi))
  sigmahat
}
# trsigma(x)

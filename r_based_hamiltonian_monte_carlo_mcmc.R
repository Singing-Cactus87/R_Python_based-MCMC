library(truncnorm)
library(coda)
library(invgamma)
library(mnormt)
library(posterior)
library(mvtnorm)
library(adaptMCMC)
#library(numDeriv)

#data generation
set.seed(321321)
dt <- rnorm(100, mean=2,sd=sqrt(0.5))

prior_mu <- function(x,a,b,log=F){
  if (log == F){
    return(dnorm(x, mean=a, sd=sqrt(b)))
  }else{
    return(dnorm(x,mean=a,sd=sqrt(b),log=T))
  }
}

prior_sig2 <- function(x,a,b,log=F){
  if (log == F){
    return(dinvgamma(x,shape=a,rate=b))
  }else{
    return(dinvgamma(x,shape=a,rate=b,log=T))
  }
}

Likelihood <- function(dt,mu,var,log=F){
  if (log==F){
    L = 1
    for (i in 1:length(dt)){
      L = L*dnorm(dt[i],mean=mu, sd=sqrt(var))
    }
    return(L)
  }else{
    L = 0
    for (i in 1:length(dt)){
      L = L+dnorm(dt[i],mean=mu, sd=sqrt(var),log=T)
    }
    return(L)
  }
}


#Leap-Frog Algorithm을 위한 derivative 정의
derivatives <- function(mu_t,sig_t,a1,b1,a2,b2,dt){
  dr1 = 0; dr2 = 0
  for (i in 1:length(dt)){
    dr1 <- dr1 +  (dt[i]-mu_t)/exp(sig_t)
    dr2 <- dr2 + (dt[i]-mu_t)^2/(2*exp(sig_t)^2)-1/(2*exp(sig_t))
  }
  dr1 <- dr1 + (a1-mu_t)/b1
  dr2 <- dr2 -(a2+1)/exp(sig_t)+b2/exp(sig_t)^2
  return(c(dr1,dr2))
}

#HMC MCMC 알고리즘
#HMC MCMC
K = 1000
t=0
mu_init = 1
sig_init = 1
momentum_init = c(0.01,0.01)
L= 10
eps = 0.004


mu_c = mu_init*1
sig_c = sig_init*1
momentum_c = momentum_init*1
mu_set = c(rep(0,K))
sig2_set = c(rep(0,K))
momentum_set1 = c(rep(0,K))
momentum_set2 = c(rep(0,K))
M = diag(2)*1
set.seed(2026)
for (i in 1:K){
  momentum_new <- rmvnorm(1,mean=c(0,0),sigma=M)
  mu_new <- mu_c*1 ; sig_new_2 <- exp(sig_c)*1; sig_new <- sig_c*1
  for (j in 1:L){
    momentum_new <- momentum_new + 0.5*eps*derivatives(mu_new,log(sig_new_2),a1=mean(dt),b1=0.1,a2=3,b2=1,dt)
    v = eps*solve(M)%*%t(momentum_new)
    mu_new = mu_new + v[1]; sig_new_2 = sig_new_2 + v[2]; sig_new <- log(sig_new_2)
    momentum_new <- momentum_new + 0.5*eps*derivatives(mu_new,log(sig_new_2),a1=mean(dt),b1=0.1,a2=3,b2=1,dt)
  }
  
  Post1 <- prior_mu(mu_new,a=mean(dt),b=0.1)*prior_sig2(exp(sig_new),a=3,b=1)*Likelihood(dt,mu=mu_new, var=exp(sig_new))
  Post2 <- prior_mu(mu_c,a=mean(dt),b=0.1)*prior_sig2(exp(sig_c),a=3,b=1)*Likelihood(dt,mu=mu_c, var=exp(sig_c))
  alpha <- (Post1*dmvnorm(c(mu_c,exp(sig_c)), mean=c(mu_new,exp(sig_new)),sigma=diag(2)*0.02)*dmvnorm(c(momentum_new),mean=c(0,0),sigma=M))/(Post2*dmvnorm(c(mu_new,exp(sig_new)), mean=c(mu_c,exp(sig_c)),sigma=diag(2)*0.02)*dmvnorm(c(momentum_c),mean=c(0,0),sigma=M))
  U <- runif(1,0,1)
  if (U < alpha){
    mu_c <- mu_new
    sig_c <- sig_new
    momentum_c <- momentum_new
    mu_set[i] <- mu_c ; sig2_set[i] <- exp(sig_c); momentum_set1[i] <- momentum_c[1]; momentum_set2[i] <- momentum_c[2]
  }else{
    mu_c <- mu_c
    sig2_c <- sig_c
    momentum_c <- momentum_c
    mu_set[i] <- mu_c ; sig2_set[i] <- exp(sig_c); momentum_set1[i] <- momentum_c[1]; momentum_set2[i] <- momentum_c[2]
  }
}

#HMC MCMC 샘플링 결과
plot.ts(mu_set)
plot.ts(sig2_set)

#burn-in 조절
n_mu_set = mu_set[1:1000]
n_sig2_set = sig2_set[1:1000]

#모수별 R_hat
print(posterior::rhat(mcmc(n_mu_set)))
print(posterior::rhat(mcmc(n_sig2_set)))


######################


#HMC_MCMC_log, 로그 사후분포, 로그 사전분포 기반 HMC 샘플링
#HMC MCMC log
K = 1000
t=0
mu_init = 1
sig_init = 1
momentum_init = c(0.01,0.01)
L= 10
eps = 0.004


mu_c = mu_init*1
sig_c = sig_init*1
momentum_c = momentum_init*1
mu_set = c(rep(0,K))
sig2_set = c(rep(0,K))
momentum_set1 = c(rep(0,K))
momentum_set2 = c(rep(0,K))
M = diag(2)*1
set.seed(2026)
for (i in 1:K){
  momentum_new <- rmvnorm(1,mean=c(0,0),sigma=M)
  mu_new <- mu_c*1 ; sig_new_2 <- exp(sig_c)*1; sig_new <- sig_c*1
  for (j in 1:L){
    momentum_new <- momentum_new + 0.5*eps*derivatives(mu_new,log(sig_new_2),a1=mean(dt),b1=0.1,a2=3,b2=1,dt)
    v = eps*solve(M)%*%t(momentum_new)
    mu_new = mu_new + v[1]; sig_new_2 = sig_new_2 + v[2]; sig_new <- log(sig_new_2)
    momentum_new <- momentum_new + 0.5*eps*derivatives(mu_new,log(sig_new_2),a1=mean(dt),b1=0.1,a2=3,b2=1,dt)
  }
  
  Post1 <- prior_mu(mu_new,a=mean(dt),b=0.1,log=T)+prior_sig2(exp(sig_new),a=3,b=1,log=T)+Likelihood(dt,mu=mu_new, var=exp(sig_new),log=T)
  Post2 <- prior_mu(mu_c,a=mean(dt),b=0.1,log=T)+prior_sig2(exp(sig_c),a=3,b=1,log=T)+Likelihood(dt,mu=mu_c, var=exp(sig_c),log=T)
  alpha <- (Post1+dmvnorm(c(mu_c,exp(sig_c)), mean=c(mu_new,exp(sig_new)),sigma=diag(2)*0.02,log=T)+dmvnorm(c(momentum_new),mean=c(0,0),sigma=M,log=T))-(Post2+dmvnorm(c(mu_new,exp(sig_new)), mean=c(mu_c,exp(sig_c)),sigma=diag(2)*0.02,log=T)+dmvnorm(c(momentum_c),mean=c(0,0),sigma=M,log=T))
  U <- runif(1,0,1)
  if (log(U) < alpha){
    mu_c <- mu_new
    sig_c <- sig_new
    momentum_c <- momentum_new
    mu_set[i] <- mu_c ; sig2_set[i] <- exp(sig_c); momentum_set1[i] <- momentum_c[1]; momentum_set2[i] <- momentum_c[2]
  }else{
    mu_c <- mu_c
    sig2_c <- sig_c
    momentum_c <- momentum_c
    mu_set[i] <- mu_c ; sig2_set[i] <- exp(sig_c); momentum_set1[i] <- momentum_c[1]; momentum_set2[i] <- momentum_c[2]
  }
}

#MCMC 샘플링 결과
plot.ts(mu_set)
plot.ts(sig2_set)

#burn_in 조절
n_mu_set = mu_set[1:1000]
n_sig2_set = sig2_set[1:1000]

#모수별 R_hat 결과
print(rhat(mcmc(n_mu_set)))
print(rhat(mcmc(n_sig2_set)))
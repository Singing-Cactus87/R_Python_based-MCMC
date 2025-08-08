install.packages("truncnorm")
install.packages("coda")
install.packages("invgamma")
install.packages("mnormt")
install.packages("posterior")
install.packages("mvtnorm")
install.packages("adaptMCMC")

library(truncnorm)
library(coda)
library(invgamma)
library(mnormt)
library(posterior)
library(mvtnorm)
library(adaptMCMC)

#data generation
set.seed(321321)
dt <- rnorm(100, mean=2,sd=sqrt(0.5))

#starting point consideration
print(mean(dt));print(var(dt))

#mu 사전분포
prior_mu <- function(x,a,b,log=F){
  if (log == F){
    return(dnorm(x, mean=a, sd=sqrt(b)))
  }else{
    return(dnorm(x,mean=a,sd=sqrt(b),log=T))
  }
}

#sigma^2 사전분포
prior_sig2 <- function(x,a,b,log=F){
  if (log == F){
    return(dinvgamma(x,shape=a,rate=b))
  }else{
    return(dinvgamma(x,shape=a,rate=b,log=T))
  }
}

#우도함수 정의(log분포 사용 여부 추가 정의)
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

#제안분포 정의
Proposal <- function(mu,sd2){
  return(rmtruncnorm(1, mean=c(mu,sd2),varcov=diag(2)*0.01,lower=c(-Inf,0)))
}


#MH algorithm
K = 1000
t=0
mu_init = 1
sig2_init = 1

mu_c = mu_init*1
sig2_c = sig2_init*1
mu_set = c(rep(0,K))
sig2_set = c(rep(0,K))

set.seed(2026)
for (i in 1:K){
  new <- Proposal(mu_c,sig2_c)
  mu_new <- new[1]; sig2_new <- new[2]
  Post1 <- prior_mu(mu_new,a=mean(dt),b=0.1)*prior_sig2(sig2_new,a=3,b=1)*Likelihood(dt,mu=mu_new, var=sig2_new)
  Post2 <- prior_mu(mu_c,a=mean(dt),b=0.1)*prior_sig2(sig2_c,a=3,b=1)*Likelihood(dt,mu=mu_c, var=sig2_c)
  alpha <- (Post1*dmtruncnorm(c(mu_c,sig2_c), mean=c(mu_new,sig2_new),varcov=diag(2)*0.01,lower=c(-Inf,0)))/(Post2*dmtruncnorm(c(mu_new,sig2_new), mean=c(mu_c,sig2_c),varcov=diag(2)*0.01,lower=c(-Inf,0)))
  U <- runif(1,0,1)
  if (U < alpha){
    mu_c <- mu_new
    sig2_c <- sig2_new
    mu_set[i] <- mu_c ; sig2_set[i] <- sig2_c
  }else{
    mu_c <- mu_c
    sig2_c <- sig2_c
    mu_set[i] <- mu_c ; sig2_set[i] <- sig2_c
  }
}

#MH MCMC 샘플링 결과
plot.ts(mu_set);plot.ts(sig2_set)
coda::HPDinterval(mcmc(sig2_set))
print(posterior::rhat(mcmc(mu_set)));print(posterior::rhat(mcmc(sig2_set)))


########################


#MH algorithm_log: 로그 사후분포, 사전분포 기반 알고리즘
K = 1000
t=0
mu_init = 1
sig2_init = 1

mu_c = mu_init*1
sig2_c = sig2_init*1
mu_set = c(rep(0,K))
sig2_set = c(rep(0,K))

set.seed(2026)
for (i in 1:K){
  new <- Proposal(mu_c,sig2_c)
  mu_new <- new[1]; sig2_new <- new[2]
  Post1 <- prior_mu(mu_new,a=mean(dt),b=0.1,log=T)+prior_sig2(sig2_new,a=3,b=1,log=T)+Likelihood(dt,mu=mu_new, var=sig2_new,log=T)
  Post2 <- prior_mu(mu_c,a=mean(dt),b=0.1,log=T)+prior_sig2(sig2_c,a=3,b=1,log=T)+Likelihood(dt,mu=mu_c, var=sig2_c,log=T)
  alpha <- (Post1+dmtruncnorm(c(mu_c,sig2_c), mean=c(mu_new,sig2_new),varcov=diag(2)*0.01,lower=c(-Inf,0),log=T))-(Post2+dmtruncnorm(c(mu_new,sig2_new), mean=c(mu_c,sig2_c),varcov=diag(2)*0.01,lower=c(-Inf,0),log=T))
  U <- runif(1,0,1)
  if (log(U) < alpha){
    mu_c <- mu_new
    sig2_c <- sig2_new
    mu_set[i] <- mu_c ; sig2_set[i] <- sig2_c
  }else{
    mu_c <- mu_c
    sig2_c <- sig2_c
    mu_set[i] <- mu_c ; sig2_set[i] <- sig2_c
  }
}

#MH MCMC 샘플링 결과
plot.ts(mu_set);plot.ts(sig2_set)
print(posterior::rhat(mcmc(mu_set)));print(posterior::rhat(mcmc(sig2_set)))
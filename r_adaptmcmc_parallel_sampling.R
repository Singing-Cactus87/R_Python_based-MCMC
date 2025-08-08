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


###adpatMCMC 적용 방법론

init.param <- c(mu=1,log_sig2=-1)

logPosterior <- function(pars){
  theta = pars
  sum(dnorm(dt,mean=theta[1],sd=exp(theta[2]),log=T))+sum(dnorm(theta[1],mean=mean(dt),sd=0.1))+sum(dnorm(exp(theta[2]),3,1))
}

samples <- MCMC(p=logPosterior,n=1e3,init=init.param,scale=c(0.01,0.01),adapt=T, acc.rate=0.35)

mu_sam <- samples$samples[,1]
sig2_sam <- exp(samples$samples[,2])

plot.ts(mu_sam)
plot.ts(sig2_sam)

print(rhat(mu_sam))
print(rhat(sig2_sam))

##########################


#Parallel Sampling(병렬 샘플링) based on adaptMCMC
samples2 <- MCMC.parallel(p=logPosterior,n=1e3,init=init.param,scale=c(0.01,0.01),adapt=T, acc.rate=0.35, n.chain = 4,n.cpu = 2)

plot.ts(samples2[[1]]$samples[,1],col="red",ylab="mu_sam",lty=1,lwd=1)
lines(samples2[[2]]$samples[,1],col="blue",lty=2,lwd=1)
lines(samples2[[3]]$samples[,1],col="green",lty=3,lwd=1)
lines(samples2[[4]]$samples[,1],col="black",lty=4,lwd=1)
legend("bottomright",legend=c("Chain 1","Chain 2","Chain 3","Chain 4"),col=c("red","blue","green","black"),lty=1:4,lwd=1)

plot.ts(exp(samples2[[1]]$samples[,2]),col="red",ylab="sig2_sam",lty=1,lwd=1)
lines(exp(samples2[[2]]$samples[,2]),col="blue",lty=2,lwd=1)
lines(exp(samples2[[3]]$samples[,2]),col="green",lty=3,lwd=1)
lines(exp(samples2[[4]]$samples[,2]),col="black",lty=4,lwd=1)
legend("bottomright",legend=c("Chain 1","Chain 2","Chain 3","Chain 4"),col=c("red","blue","green","black"),lty=1:4,lwd=1)
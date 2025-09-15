rm(list=ls())

source("R/simulateCSCP.R")
source("R/fitCSCP_lm.R")

# barlambda <- 1000          # desired E[N]/|W|
# w         <- 1          # fraction of mean carried by stochastic CSCP part
# alphas    <- 1#c(0.5, 0.5)  # split across components (sums to 1)
# scales    <- 0.08#c(0.02, 0.20)
# 
# temp <- sample(1:1000,1)
# set.seed(temp)
# X <- simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = 1)
# plot(X)
# npoints(X)
# 
# out <- fitCSCP(X,plot=T,correction="Ripley",zerocor="JonesFoster",divisor="a",bw=bw.pcf(X,cv.method="leastSQ"))
# out$par


set.seed(123)
nsim <- 100
scales <- 0.1
alphas <- 1
w <- 1
barlambda <- 1000
Xmany <- simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = nsim)

t1 <- Sys.time()
pb <- txtProgressBar(0,nsim)
result <- matrix(NA,nsim,2)
for(i in 1:nsim){
  out <- fitCSCP(Xmany[[i]],plot=T,correction="Ripley",zerocor="JonesFoster",divisor="a",bw=bw.pcf(X,cv.method="leastSQ"))
  result[i,] <- out$par
  setTxtProgressBar(pb,i)
}
close(pb)
t2 <- Sys.time()
t2-t1

hist(result[,1],main="phi")
abline(v=scales[1],col=2,lty=2)



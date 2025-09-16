rm(list=ls())

source("R/simulateCSCP.R")
source("R/fitCSCP_lm.R")

barlambda <- 1000          # desired E[N]/|W|
w         <- 1          # fraction of mean carried by stochastic CSCP part
alphas    <- 1#c(0.5, 0.5)  # split across components (sums to 1)
scales    <- 0.08#c(0.02, 0.20)

temp <- 228#sample(1:1000,1)
set.seed(temp)
X <- simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = 1)
plot(X)
npoints(X)

out1 <- fitCSCP(X,plot=TRUE,
                correction="Ripley",zerocor="JonesFoster",divisor="a",bw=bw.pcf(X,cv.method="leastSQ"))
out1$par

out2 <- fitCSCP(X,plot=TRUE,fitalpha=FALSE,
                correction="Ripley",zerocor="JonesFoster",divisor="a",bw=bw.pcf(X,cv.method="leastSQ"))
out2$par



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
  out <- fitCSCP(Xmany[[i]],plot=F,correction="Ripley",zerocor="JonesFoster",divisor="a",bw=bw.pcf(X,cv.method="leastSQ"))
  result[i,] <- out$par
  setTxtProgressBar(pb,i)
}
close(pb)
t2 <- Sys.time()
t2-t1

hist(result[,1],main="phi")
abline(v=scales[1],col=2,lty=2)

hist(result[,2],main="(derived) alpha")
abline(v=alphas[1],col=2,lty=2)



set.seed(123)
nsim <- 100
scales <- 0.1
alphas <- 1
w <- 1
barlambda <- 1000
Xmany <- simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = nsim)

t1 <- Sys.time()
pb <- txtProgressBar(0,nsim)
result2 <- matrix(NA,nsim,2)
for(i in 1:nsim){
  out <- fitCSCP(Xmany[[i]],plot=F,fitalpha=F,correction="Ripley",zerocor="JonesFoster",divisor="a",bw=bw.pcf(X,cv.method="leastSQ"))
  result2[i,] <- out$par
  setTxtProgressBar(pb,i)
}
close(pb)
t2 <- Sys.time()
t2-t1

hist(result2[,1],main="phi")
abline(v=scales[1],col=2,lty=2)


r1mse <- sqrt(mean((result[,1]-scales[1])^2))
pdf("phi_hist1.pdf",10,5)
par(mfrow=c(1,2))
hist(result[,1],main=paste0("phi (free intercept); RMSE=",round(r1mse,3)),xlab="")
abline(v=scales[1],col=2,lty=2,lwd=3)
hist(result[,2],main="alpha (derived)",xlab="")
abline(v=alphas[1],col=2,lty=2,lwd=3)
dev.off()

r2mse <- sqrt(mean((result2[,1]-scales[1])^2))
pdf("phi_hist2.pdf",5,5)
hist(result2[,1],main=paste0("phi (locked intercept); RMSE=",round(r2mse,3)),xlab="")
abline(v=scales[1],col=2,lty=2,lwd=3)
dev.off()

source("R/simulateCSCP.R")
plot(simulateCSCP(list(list(mu=0, scale=0.02, var=500),list(mu=0, scale=0.2, var=500)), nsim=6),
     main="k = 2; mu = 0,0; scale = 0.02,0.2; var = 500,500")



## w = 1 

barlambda <- 100          # desired E[N]/|W|
w         <- 1          # fraction of mean carried by stochastic CSCP part
alphas    <- c(0.5, 0.5)  # split across components (sums to 1)
scales    <- c(0.02, 0.20)

par(mfrow = c(2,3))
plot(simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = 6),
     main = expression(paste("reparam CSCP:  ", bar(lambda),"=100;  w=1;  ",alpha,"=0.5,0.5;  ",phi,"=0.02,0.2")))

barlambda <- 10000          # desired E[N]/|W|
w         <- 1          # fraction of mean carried by stochastic CSCP part
alphas    <- c(0.5, 0.5)  # split across components (sums to 1)
scales    <- c(0.02, 0.20)
par(mfrow = c(2,3))
plot(simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = 6),
     main = expression(paste("reparam CSCP:  ", bar(lambda),"=10000;  w=1;  ",alpha,"=0.5,0.5;  ",phi,"=0.02,0.2")))


## w < 1

barlambda <- 1000          # desired E[N]/|W|
w         <- 0.75         # fraction of mean carried by stochastic CSCP part
alphas    <- c(0.5, 0.5)  # split across components (sums to 1)
scales    <- c(0.02, 0.20)

par(mfrow = c(2,3))
plot(simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = 6),
     main = expression(paste("reparam CSCP:  ", bar(lambda),"=1000;  w=0.75;  ",alpha,"=0.5,0.5;  ",phi,"=0.02,0.2")))

barlambda <- 1000          # desired E[N]/|W|
w         <- 0.1          # fraction of mean carried by stochastic CSCP part
alphas    <- c(0.5, 0.5)  # split across components (sums to 1)
scales    <- c(0.02, 0.20)
par(mfrow = c(2,3))
plot(simulateCSCP_reparam(barlambda, w, alphas, scales, nsim = 6),
     main = expression(paste("reparam CSCP:  ", bar(lambda),"=1000;  w=0.1;  ",alpha,"=0.5,0.5;  ",phi,"=0.02,0.2")))


library(fda)
### setting ###
z <- 1:100/100; Rep <- 200
N1 <- 1000; N2 <- 700; 
N.test <- 100;
N1.train <- N1-N.test
N2.train <- N2-N.test
basis <- create.bspline.basis(c(0,1),nbasis=24)
a <- sqrt(40); b <- sqrt((100-a^2)/7)
#### simulate functions ####
# for B-spline #
sigma1 <- rep(c(b,b,a,rep(b,5)),3)
sigma2 <- rep(c(rep(b,5),a,b,b),3)

p.vpc1 <- vector(length=Rep)
p.vpc2 <- vector(length=Rep)
for(rep in 1:Rep){
  score1 <- matrix(rnorm(N1*24,0,sigma1),24,N1)
  score2 <- matrix(rnorm(N2*24,0,sigma2),24,N2)
  F1 <- fd(score1,basis)
  F2 <- fd(score2,basis)
  F1.dis <- eval.fd(z,F1)
  F2.dis <- eval.fd(z,F2)
  C1 <- cov(t(F1.dis[,1:N1.train]))
  C2 <- cov(t(F2.dis[,1:N2.train]))
  ### VPC ###
  cat.vpc1 <- vector(length=100)
  cat.vpc2 <- vector(length=100)
  C.diff <- (C1-C2)%*%(C1-C2)
  eigen.result <- eigen(C.diff)
  d <- which(cumsum(eigen.result$values)/sum(eigen.result$values)>0.9)[1]
  ff <- eigen.result$vector[,1:d]*sqrt(length(z)) # feature function
  s1 <- t(F1.dis[,-(1:N1.train)])%*%ff/length(z)
  s2 <- t(F2.dis[,-(1:N2.train)])%*%ff/length(z)
  S1 <- t(ff)%*%C1%*%ff/length(z)^2
  S2 <- t(ff)%*%C2%*%ff/length(z)^2
  for(i in 1:100){
    diff1 <- sum((matrix(s1[i,],ncol=1)%*%s1[i,]-S1)^2)
    diff2 <- sum((matrix(s1[i,],ncol=1)%*%s1[i,]-S2)^2)
    if(diff1<diff2) cat.vpc1[i] <- 1; if(diff1>diff2) cat.vpc1[i] <- 2
    diff1 <- sum((matrix(s2[i,],ncol=1)%*%s2[i,]-S1)^2)
    diff2 <- sum((matrix(s2[i,],ncol=1)%*%s2[i,]-S2)^2)
    if(diff1<diff2) cat.vpc2[i] <- 1; if(diff1>diff2) cat.vpc2[i] <- 2
  }
  p.vpc1[rep] <- sum(cat.vpc1==1)/100
  p.vpc2[rep] <- sum(cat.vpc2==2)/100
}

print(c(mean(p.vpc1),sd(p.vpc1),mean(p.vpc2),sd(p.vpc2)))

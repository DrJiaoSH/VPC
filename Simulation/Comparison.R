library(fda)
### setting ###
z <- 1:100/100; Rep <- 200
N1 <- 200; N2 <- 200; N.train <- 100
basis <- create.bspline.basis(c(0,1),nbasis=24)
a <- sqrt(80); b <- sqrt((100-a^2)/7)
#### simulate functions ####
sigma1 <- rep(c(b,b,a,rep(b,5)),3)
sigma2 <- rep(c(rep(b,5),a,b,b),3)

p.vpc1 <- vector(length=Rep)
p.vpc2 <- vector(length=Rep)
p.pj1 <- vector(length=Rep)
p.pj2 <- vector(length=Rep)
p.fqc1 <- vector(length=Rep)
p.fqc2 <- vector(length=Rep)
for(rep in 1:Rep){
  print(rep)
  score1 <- matrix(rnorm(N1*24,0,sigma1),24,N1)
  score1 <- score1-rowMeans(score1)
  score2 <- matrix(rnorm(N2*24,0,sigma2),24,N2)
  score2 <- score2-rowMeans(score2)
  F1 <- fd(score1,basis)
  F2 <- fd(score2,basis)
  F1.dis <- eval.fd(z,F1)
  F2.dis <- eval.fd(z,F2)
  C1 <- cov(t(F1.dis[,1:N.train]))
  C2 <- cov(t(F2.dis[,1:N.train]))
 
  # This section is for the VPC mathod
  cat.vpc1 <- vector(length=100)
  cat.vpc2 <- vector(length=100)
  C.diff <- (C1-C2)%*%(C1-C2)
  eigen.result <- eigen(C.diff)
  d <- which(cumsum(eigen.result$values)/sum(eigen.result$values)>0.9)[1]
  ff <- eigen.result$vector[,1:d]*sqrt(length(z)) # feature function
  s1 <- t(F1.dis[,-(1:N.train)])%*%ff/length(z)
  s2 <- t(F2.dis[,-(1:N.train)])%*%ff/length(z)
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
  
  # This method is for the projection method
  eigen.result1 <- eigen(C1)
  eigen.result2 <- eigen(C2)
  d1 <- which(cumsum(eigen.result1$values)/sum(eigen.result1$values)>0.9)[1]
  d2 <- which(cumsum(eigen.result2$values)/sum(eigen.result2$values)>0.9)[1]
  fpc1 <- eigen(C1)$vector[,1:d1]*sqrt(length(z))
  fpc2 <- eigen(C2)$vector[,1:d2]*sqrt(length(z))
  proj1 <- t(F1.dis[,-(1:N.train)])%*%fpc1%*%t(fpc1)/length(z)
  proj2 <- t(F1.dis[,-(1:N.train)])%*%fpc2%*%t(fpc2)/length(z)
  diff1 <- apply((F1.dis[,-(1:N.train)] - t(proj1))^2,2,mean)
  diff2 <- apply((F1.dis[,-(1:N.train)] - t(proj2))^2,2,mean)
  cat.pj1 <- sum(diff1<diff2)
  proj1 <- t(F2.dis[,-(1:N.train)])%*%fpc1%*%t(fpc1)/length(z)
  proj2 <- t(F2.dis[,-(1:N.train)])%*%fpc2%*%t(fpc2)/length(z)
  diff1 <- apply((F2.dis[,-(1:N.train)] - t(proj1))^2,2,mean)
  diff2 <- apply((F2.dis[,-(1:N.train)] - t(proj2))^2,2,mean)
  cat.pj2 <- sum(diff1>diff2)
  p.pj1[rep] <- cat.pj1/100
  p.pj2[rep] <- cat.pj2/100
  
  # This is for the likehood ratio test
  cat.fqc1 <- vector(length=100)
  cat.fqc2 <- vector(length=100)
  eigent <- eigen(cov(t(cbind(F1.dis[,1:N.train],F2.dis[,1:N.train]))))
  d.fqc <- which(cumsum(eigent$values)/sum(eigent$values)>0.9)[1]
  eigenfun <- eigent$vectors[,1:d.fqc]*sqrt(length(z))
  eigenv1 <- diag(t(eigenfun)%*%C1%*%eigenfun)/length(z)^2
  eigenv2 <- diag(t(eigenfun)%*%C2%*%eigenfun)/length(z)^2
  escore1 <- t(F1.dis[,-(1:N.train)])%*%eigenfun/length(z)
  escore2 <- t(F2.dis[,-(1:N.train)])%*%eigenfun/length(z)
  for(i in 1:100){
    Q <- 0.5*(sum(log(eigenv1)-log(eigenv2))-sum(escore1[i,]^2/eigenv2-escore1[i,]^2/eigenv1))
    if(Q>0) cat.fqc1[i] <- 2; if(Q<0) cat.fqc1[i] <- 1
    Q <- 0.5*(sum(log(eigenv1)-log(eigenv2))-sum(escore2[i,]^2/eigenv2-escore2[i,]^2/eigenv1))
    if(Q>0) cat.fqc2[i] <- 2; if(Q<0) cat.fqc2[i] <- 1
  }
  p.fqc1[rep] <- sum(cat.fqc1==1)/100
  p.fqc2[rep] <- sum(cat.fqc2==2)/100
}


library(fda)
### setting ###
z <- 1:100/100; Rep <- 200
N1 <- 300; N2 <- 300; N.train <- 200
basis <- create.fourier.basis(c(0,1),nbasis=21)

# for B-spline #
sigma1 <- c(1,1,1,0,0,rep(1,16))
sigma2 <- c(1,0,0,1,1,rep(1,16))

sigma1 <- c(1,1,0,1,0,1,0,1,rep(1,13))
sigma2 <- c(1,0,1,0,1,0,1,1,rep(1,13))

p.vpc1 <- list(); p.vpc2 <- list()
p.vpc1.inf <- vector(length=Rep)
p.vpc2.inf <- vector(length=Rep)
for(rep in 1:Rep){
  p.vpc1[[Rep]] <- vector(length=9)
  p.vpc2[[Rep]] <- vector(length=9)
  score1 <- matrix(rnorm(N1*21,0,sigma1),21,N1)
  score2 <- matrix(rnorm(N2*21,0,sigma2),21,N2)
  F1 <- fd(score1,basis)
  F2 <- fd(score2,basis)
  F1.dis <- eval.fd(z,F1)
  F2.dis <- eval.fd(z,F2)
  C1 <- cov(t(F1.dis[,1:N.train]))
  C2 <- cov(t(F2.dis[,1:N.train]))
  C.diff <- (C1-C2)%*%(C1-C2)
  eigen.result <- eigen(C.diff)
  for(d in 1:9){
    cat.vpc1 <- vector(length=100)
    cat.vpc2 <- vector(length=100)
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
    p.vpc1[[rep]][d] <- sum(cat.vpc1==1)/100
    p.vpc2[[rep]][d] <- sum(cat.vpc2==2)/100
  }
  
  cat.vpc1 <- vector(length=100)
  cat.vpc2 <- vector(length=100)
  for(i in 1:100){
    diff1 <- sum((C1-F1.dis[,N.train+i]%*%t(F1.dis[,N.train+i]))^2)
    diff2 <- sum((C2-F1.dis[,N.train+i]%*%t(F1.dis[,N.train+i]))^2)
    if(diff1<diff2) cat.vpc1[i] <- 1; if(diff1>diff2) cat.vpc1[i] <- 2
    diff1 <- sum((C1-F2.dis[,N.train+i]%*%t(F2.dis[,N.train+i]))^2)
    diff2 <- sum((C2-F2.dis[,N.train+i]%*%t(F2.dis[,N.train+i]))^2)
    if(diff1<diff2) cat.vpc2[i] <- 1; if(diff1>diff2) cat.vpc2[i] <- 2
  }
  p.vpc1.inf[rep] <- sum(cat.vpc1==1)/100
  p.vpc2.inf[rep] <- sum(cat.vpc2==2)/100
}

#P.vpc1.case1 <- c(colMeans(do.call(rbind,p.vpc1)),mean(p.vpc1.inf))
#P.vpc2.case1 <- c(colMeans(do.call(rbind,p.vpc2)),mean(p.vpc2.inf))

#P.vpc1.case2 <- c(colMeans(do.call(rbind,p.vpc1)),mean(p.vpc1.inf))
#P.vpc2.case2 <- c(colMeans(do.call(rbind,p.vpc2)),mean(p.vpc2.inf))

## ff ##
sigma1 <- c(1,1,1,0,0,rep(1,16))
sigma2 <- c(1,0,0,1,1,rep(1,16))

sigma1 <- c(1,1,0,1,0,1,0,1,rep(1,13))
sigma2 <- c(1,0,1,0,1,0,1,1,rep(1,13))

score1 <- matrix(rnorm(N1*21,0,sigma1),21,N1)
score2 <- matrix(rnorm(N2*21,0,sigma2),21,N2)
F1 <- fd(score1,basis)
F2 <- fd(score2,basis)
F1.dis <- eval.fd(z,F1)
F2.dis <- eval.fd(z,F2)
C1 <- cov(t(F1.dis))
C2 <- cov(t(F2.dis))
C.diff <- (C1-C2)%*%(C1-C2)
eigen.result <- eigen(C.diff)
ff <- eigen.result$vectors[,1:9]

#ff.case1 <- ff
#ff.case2 <- ff
## plot ##
par(mfrow=c(2,2),mar=c(2,2,2,2))

plot(P.vpc1.case1,ylim=c(0.5,1),xaxt="n",ylab="Classification rate (%)",xlab="Dimension",main="ACR (Setting 1)")
axis(1,at=1:10,label=c(1:9,"Inf"))
lines(P.vpc2.case1,type="p",pch=3)
legend("bottomright",legend=c("Group 0", "Group 1"),pch = c(1,3))

plot(ff.case1[,1],xaxt="n",ylim=c(-0.25,0.25),type="l",main="Discrepancy feature functions",ylab="",xlab="Time")
for(i in 2:4) lines(ff.case1[,i],lty=i)
axis(1,at=c(1,20,40,60,80,100),label=c(0,0.2,0.4,0.6,0.8,1.0))

plot(P.vpc1.case2,ylim=c(0.5,1),xaxt="n",ylab="Classification rate (%)",xlab="Dimension",main="ACR (Setting 2)")
axis(1,at=1:10,label=c(1:9,"Inf"))
lines(P.vpc2.case2,type="p",pch=3)
legend("bottomright",legend=c("Group 0", "Group 1"),pch = c(1,3))

plot(ff.case2[,1],xaxt="n",ylim=c(-0.25,0.25),type="l",main="Discrepancy feature functions",ylab="",xlab="Time")
for(i in 2:6) lines(ff.case2[,i],lty=i)
axis(1,at=c(1,20,40,60,80,100),label=c(0,0.2,0.4,0.6,0.8,1.0))

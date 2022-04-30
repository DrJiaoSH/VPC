library(fda)

##### simulate FMA processes #####
rarmah22 = function(n, basis, Phi1, Phi2, Phi3, Phi4, Phi5, Phi6, Theta1, Theta2, Sigma){
  # build coefficients matrix 
  noise = matrix(rnorm(D*n,0,Sigma),D,n)
  coef = noise
  
  # recursion
  for (i in 7:n)
    coef[,i] = Phi1 %*% noise[,i-1] + Phi2 %*% noise[,i-2] +
               Phi3 %*% noise[,i-3] + Phi4 %*% noise[,i-4] +
               Phi5 %*% noise[,i-5] + Phi6 %*% noise[,i-6] + noise[,i]
  
  fd(coef, basis=basis)
}

D = 21; n <- 600 # n is the sample size, D is the number of Fourier basis

# Specify the covariance structure of the two groups
Sigma1 = c(1,rep(c(0.8,0.8,1,1),5))
Sigma2 = c(1,rep(c(1,1,0.8,0.8),5))
four = create.fourier.basis(rangeval=c(0, 1), nbasis=D)
Psi1 = matrix(0,D,D); Psi2 = matrix(0,D,D)
for(i in 1:D){for(j in 1:D){Psi1[i,j]=rnorm(1,0,Sigma1[i]*Sigma1[j])}}
for(i in 1:D){for(j in 1:D){Psi2[i,j]=rnorm(1,0,Sigma2[i]*Sigma2[j])}}
Psi1 = Psi1/Mod(eigen(Psi1)$values[1])
Psi2 = Psi2/Mod(eigen(Psi2)$values[1])

# Specify the coefficient operator of the FMA process
a1 <- 0.4; a2 <- 0.4
Phi11 = a1*Psi1; Phi12 = a2*Psi2
a3 <- 0.4; a4 <- 0.4
Phi21 = a3*Psi1; Phi22 = a4*Psi2
a5 <- 0.4; a6 <- 0.4
Phi31 = a5*Psi1; Phi32 = a6*Psi2
a7 <- 0.0; a8 <- 0.0
Phi41 = a7*Psi1; Phi42 = a8*Psi2
a9 <- 0.0; a10 <- 0.0
Phi51 = a9*Psi1; Phi52 = a10*Psi2
a11 <- 0.0; a12 <- 0.0
Phi61 = a11*Psi1; Phi62 = a12*Psi2

z <- 1:100/100; p <- 5; ny <- 100; Rep <- 200 # p is the maximum lag, z is the grids, Rep is the repetition

index <- list()
for(k in 1:ny){
  index[[k]] <- (p*(k-1)+1):(p*k)
}
prop1 <- matrix(0,Rep,p); prop2 <- matrix(0,Rep,p)
s <- 10
#### classification #####
for(rep in 1:Rep){
  ### two groups of functions, FTS1 and FTS2 ###
  FTS1=rarmah22(n=n, basis=four, Phi1=Phi11, Phi2=Phi21, Phi3=Phi31, Phi4=Phi41, Phi5=Phi51, Phi6=Phi61, Sigma=Sigma1)
  FTS2=rarmah22(n=n, basis=four, Phi1=Phi12, Phi2=Phi22, Phi3=Phi32, Phi4=Phi42, Phi5=Phi52, Phi6=Phi62, Sigma=Sigma2)

  ### scale the functions to achieve unbalanced classification ###
  fts1 <- eval.fd(z,FTS1); fts2 <- eval.fd(z,FTS2)
  fts1 <- fts1 %*% diag(1/sqrt(diag(t(fts1)%*%fts1/length(z))))
  fts2 <- fts2 %*% diag(1/sqrt(diag(t(fts2)%*%fts2/length(z))))
  
  ### calculate the (lagged) second moment
  G1 <- list(); G2 <- list()
  for(i in 1:p){
    G1[[i]] <- fts1[,1:(n-i+1)]%*%t(fts1[,i:n])/(n-i+1)
    G1[[i]] <- G1[[i]]+t(G1[[i]])
    G2[[i]] <- fts2[,1:(n-i+1)]%*%t(fts2[,i:n])/(n-i+1)
    G2[[i]] <- G2[[i]]+t(G2[[i]])
  }
  
  ### norm of different lags ###
  W <- vector(length=p) 
  for(i in 1:p){
    W[i] <- ((sqrt(sum(G1[[i]]^2))+sqrt(sum(G2[[i]]^2)))/2)^(-1)
  }
  
  ### obtain the discriminative feature functions (df) and dimension (d) ###
  df <- list(); d <- vector(length=p) 
  for(i in 1:p){
    pca <- eigen((G1[[i]]- G2[[i]])%*%(G1[[i]]- G2[[i]]))
    d[i] <- which(cumsum(pca$values)/sum(pca$values)>0.9)[1]
    df[[i]] <- eigen((G1[[i]]- G2[[i]])%*%(G1[[i]]- G2[[i]]))$vectors[,1:d[i]]
  }
  
  ### obtain the scores ###
  S1 <- list(); S2 <- list() 
  for(i in 1:p){
    S1[[i]] <- t(df[[i]])%*%G1[[i]]%*%df[[i]]/length(z)^2
    S2[[i]] <- t(df[[i]])%*%G2[[i]]%*%df[[i]]/length(z)^2
  }
  
  ##### New data generation ######
  Y1 <- rarmah22(n=ny*p, basis=four, Phi1=Phi11, Phi2=Phi21, Phi3=Phi31, Phi4=Phi41, Phi5=Phi51, Phi6=Phi61, Sigma=Sigma1)
  Y2 <- rarmah22(n=ny*p, basis=four, Phi1=Phi12, Phi2=Phi22, Phi3=Phi32, Phi4=Phi42, Phi5=Phi52, Phi6=Phi62, Sigma=Sigma2)
  y1 <- eval.fd(z,Y1)
  y2 <- eval.fd(z,Y2)
  S1.y <- list()
  for(k in 1:ny){
    S1.y[[k]] <- list()
    for(i in 1:p){  
      S1.y[[k]][[i]] <- t(df[[i]])%*%(y1[,index[[k]][1:(p-i+1)]]%*%t(y1[,index[[k]][i:p]])/(p-i+1))%*%df[[i]]/length(z)^2
                       +t(df[[i]])%*%(y1[,index[[k]][i:p]]%*%t(y1[,index[[k]][1:(p-i+1)]])/(p-i+1))%*%df[[i]]/length(z)^2
    }
  }
  S2.y <- list()
  for(k in 1:ny){
    S2.y[[k]] <- list()
    for(i in 1:p){  
      S2.y[[k]][[i]] <- t(df[[i]])%*%(y2[,index[[k]][1:(p-i+1)]]%*%t(y2[,index[[k]][i:p]])/(p-i+1))%*%df[[i]]/length(z)^2
                       +t(df[[i]])%*%(y2[,index[[k]][i:p]]%*%t(y2[,index[[k]][1:(p-i+1)]])/(p-i+1))%*%df[[i]]/length(z)^2
    }
  }
  
  ##### calculate Hilbert-Schimidt metric #####
  dist11 <- list(); dist12 <- list()
  for(k in 1:ny){
    dist11[[k]] <- vector(length=p); dist12[[k]] <- vector(length=p)
    for(i in 1:p){
      dist11[[k]][i] <- sum((S1.y[[k]][[i]]-S1[[i]])^2)
      dist12[[k]][i] <- sum((S1.y[[k]][[i]]-S2[[i]])^2)
    }
  }
  dist21 <- list(); dist22 <- list()
  for(k in 1:ny){
    dist21[[k]] <- vector(length=p); dist22[[k]] <- vector(length=p)
    for(i in 1:p){
      dist21[[k]][i] <- sum((S2.y[[k]][[i]]-S1[[i]])^2)
      dist22[[k]][i] <- sum((S2.y[[k]][[i]]-S2[[i]])^2)
    }
  }

  ##### classification rate of single lag #####
  di1 <- list(); di2 <- list(); w1 <- vector(length=p)
  for(l in 1:p){
    di1[[l]] <- vector(length=ny); di2[[l]] <- vector(length=ny)
    for(k in 1:ny){
      di1[[l]][k] <- sum(W[l]*dist11[[k]][l])
      di2[[l]][k] <- sum(W[l]*dist12[[k]][l])
    }
    w1[l] <- sum(di1[[l]]<di2[[l]])/ny
  }
  di1 <- list(); di2 <- list(); w2 <- vector(length=p)
  for(l in 1:p){
    di1[[l]] <- vector(length=ny); di2[[l]] <- vector(length=ny)
    for(k in 1:ny){
      di1[[l]][k] <- sum(W[l]*dist21[[k]][l])
      di2[[l]][k] <- sum(W[l]*dist22[[k]][l])
    }
    w2[l] <- sum(di1[[l]]>di2[[l]])/ny
  }
  w <- (w1+w2)/2 # classification rate based on different lags

  #### classification ####
  W.new <- W*exp(s*w) # weight function W(h)
  Dist1 <- vector(length=ny); Dist2 <- vector(length=ny)
  for(l in 1:p){
    for(k in 1:ny){
      Dist1[k] <- sum(W.new[1:l]*dist11[[k]][1:l])
      Dist2[k] <- sum(W.new[1:l]*dist12[[k]][1:l])
    }
    prop1[rep,l] <- sum(Dist1<Dist2)/ny # classification rate of group 1
  }
  Dist1 <- vector(length=ny); Dist2 <- vector(length=ny)
  for(l in 1:p){
    for(k in 1:ny){
      Dist1[k] <- sum(W.new[1:l]*dist21[[k]][1:l])
      Dist2[k] <- sum(W.new[1:l]*dist22[[k]][1:l])
    }
    prop2[rep,l] <- sum(Dist1>Dist2)/ny # classification rate of group 2
  }
}

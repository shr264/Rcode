library(caret)
library(lassoshooting)
library(Matrix)
library(MASS)
library(lassoshooting)
library(microbenchmark)
library(emulator)
library(ROCR)
library(gridExtra)
library(reshape)
library(ggplot2)
library(caret)
library(xtable)
source('CSCS2.r')
source('huang.r')
source('shoaje.r')

load("raw.RData")
X.test = scale(X.test, center = FALSE, scale = TRUE)
X.train = scale(X.train, center = FALSE, scale = TRUE)

BIC <- function(n,Lhat,S){
    Ohat = t(Lhat)%*%Lhat
    bic = n*(-log(det(Ohat))+sum(diag(Ohat%*%S))) + log(n)*sum(Lhat[lower.tri(Lhat)]!=0)
    return(bic)
}

methodfromname <- function(X, lambda, name){
    if(name == "cscs"){
        out = CSCS2(X, lambda)}
    else if(name == "huang"){
        out = huang(X, lambda, returnL=TRUE)}
    else if(name == "shoj"){
        out = shoaje(X, lambda)
    }
    else{ print('Method is incorrect')
          out = NULL}
    return(out)
}

outputfromname <- function(out1, name){
    if(name == "cscs"){
        out = out1$L}
    else if(name == "huang"){
        out = out1$L}
    else if(name == "shoj"){
        out = out1$T
    }
    else{print('Method is incorrect')
          out = NULL}
    return(out)
}

cv = function(lambda,K,name,data){
    cvec = rep(0,K)
    set.seed(12345)
    foldsize <- floor(dim(data)[1]/K)
    rndmsmpl <- sample(1:dim(data)[1])
    folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize))
    for(k in 1:K){
        Ytrain = data[-c(folds[[k]]),]
        Ytest = data[c(folds[[k]]),]
        E = scale(Ytrain, center = TRUE, scale = FALSE)
        out1 = methodfromname(as.matrix(E), lambda, name)
        out = outputfromname(out1,name)
        Sigma_v = t(out)%*%out
        s_v = dim(Ytest)[1]
        sumterm = sum(diag(Ytest%*%(Sigma_v)%*%t(Ytest)))
        cvec[k] = s_v*log(det(solve(Sigma_v)))+sumterm
    }
    cv = (1/K)*sum(cvec)
    return(cv)
}

loglik <- function(L,data,mu){
    S <- cov(data)
    Xbar <- colMeans(data)
    n <- nrow(data)
    Omega <- t(L)%*%L
    Sigma <- solve(Omega)
    out <- -(n/2)*(log(det(Sigma)) + sum(diag(Omega%*%S)) + t(Xbar-mu)%*%Omega%*%(Xbar-mu))
    return(out)
}

mu <- colMeans(X.train)

name = "cscs"
lambdal = 40
lambda = seq(0.01,0.4,length=lambdal)
out = rep(0,length(lambda))
for(k in 1:length(lambda)){
    out[k] = cv(lambda[k], 5 ,name, X.train)
}
lambdastar = lambda[which.min(out)]
E = X.train
out1 = methodfromname(as.matrix(E), lambdastar, name)
out = outputfromname(out1,name)
Shat = solve(t(out)%*%out)
CSCS = Shat
cscsll <- loglik(out,X.test,mu)

name = "huang"
lambda = seq(0.01,0.8,length=lambdal)
out = rep(0,length(lambda))
for(k in 1:length(lambda)){
    out[k] = cv(lambda[k], 5 ,name, X.train)
}
lambdastar = lambda[which.min(out)]
E = X.train
out1 = methodfromname(as.matrix(E), lambdastar, name)
out = outputfromname(out1,name)
Shat = solve(t(out)%*%out)
huang = Shat
huangll <- loglik(out,X.test,mu)

name = "shoj"
lambda = seq(0.01,0.4,length=lambdal)
out = rep(0,length(lambda))
for(k in 1:length(lambda)){
    out[k] = cv(lambda[k], 5 ,name, X.train)
}
lambdastar = lambda[which.min(out)]
E = X.train
out1 = methodfromname(as.matrix(E), lambdastar, name)
out = outputfromname(out1,name)
Shat = solve(t(out)%*%out)
shoj = Shat
shojll <- loglik(out,X.test,mu)

ll <- c(cscsll,huangll,shojll)
names(ll) <- c("cscs","huang","shoj")

write(ll,file = paste("ll.txt",sep=""))

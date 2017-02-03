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
source('CSCS2.r')
source('huang.r')
source('shoaje.r')

loglik <- function(L,data,mu){
    S <- cov(data)
    Xbar <- colMeans(data)
    n <- nrow(data)
    Omega <- t(L)%*%L
    Sigma <- solve(Omega)
    out <- -(n/2)*(log(det(Sigma)) + sum(diag(Omega%*%S)) + t(Xbar-mu)%*%Omega%*%(Xbar-mu))
    return(out)
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

trainsize = 205
ccd = as.matrix(read.table("call-center-data.txt",header=FALSE))
dim(ccd)
### in pourhamadi et al, they apply the transformation sqrt(N_it + 1/4). Has this already been done? No! We apply the transformation here.
data = sqrt(ccd + 1/4)
data = scale(data,center=TRUE,scale=FALSE)
set.seed(12345)
trainindex = 1:trainsize; ord = 'ordered'
###trainindex = sample(1:nrow(data), size = trainsize); ord = 'random'
train = data[trainindex,]
mu <- colMeans(train)
test = data[-trainindex,]
time = 52:102

cv = function(lambda,K,name){
    cvec = rep(0,K)
    set.seed(12345)
    index = sample(1:trainsize,replace=FALSE)
    foldsize = ceiling(trainsize/K)
    inTrain <- list()
    for(i in 1:(K-1)){
       inTrain[[i]] <- index[((i-1)*foldsize + 1): ((i-1)*foldsize + foldsize)]
   }
    inTrain[[K]] <- index[((K-1)*foldsize + 1):length(index)]
    for(k in 1:K){
        Ytrain = train[-c(inTrain[[k]]),]
        Ytest = train[c(inTrain[[k]]),]
        E = Ytrain
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

Shat = (t(train)%*%train)/dim(train)[2]
S = Shat

if(trainsize>102){
    R = chol(S)
    L = t(solve(R))
    sll <- loglik(L,test,mu)
}

name = "cscs"
lambdal = 0.01
lambda = seq(0.01,0.4,by=lambdal)
out = rep(0,length(lambda))
for(k in 1:length(lambda)){
    out[k] = cv(lambda[k], 5 ,name)
}
plot(out)
which.min(out)
lambdastar = lambda[which.min(out)]
E = train
out1 = methodfromname(as.matrix(E), lambdastar, name)
out = outputfromname(out1,name)
Shat = solve(t(out)%*%out)
CSCS = Shat
cscsll <- loglik(out,test,mu)


name = "huang"
lambda = seq(0.3,0.8,by = lambdal)
out = rep(0,length(lambda))
for(k in 1:length(lambda)){
    out[k] = cv(lambda[k], 5 ,name)
}
plot(out)
which.min(out)
lambdastar = lambda[which.min(out)]
E = train
out1 = methodfromname(as.matrix(E), lambdastar, name)
out = outputfromname(out1,name)
Shat = solve(t(out)%*%out)
Huang = Shat
huangll <- loglik(out,test,mu)

name = "shoj"
lambda = seq(0.01,0.4,by = lambdal)
out = rep(0,length(lambda))
for(k in 1:length(lambda)){
    out[k] = cv(lambda[k], 5 ,name)
}
plot(out)
which.min(out)
lambdastar = lambda[which.min(out)]
E = train
out1 = methodfromname(as.matrix(E), lambdastar, name)
out = outputfromname(out1,name)
Shat = solve(t(out)%*%out)
Shoj = Shat
shojll <- loglik(out,test,mu)


predict.mean <- function(x1,mu,Sigma){
    p1 <- length(x1)
    p <- length(mu)
    p2 <- p-p1
    mu1 <- mu[1:p1]
    mu2 <- mu[(p1+1):p]
    Sigma11 <- Sigma[1:p1,1:p1]
    #Sigma12 <- Sigma[1:p1,(p1+1):p]
    Sigma21 <- Sigma[(p1+1):p,1:p1]
    #Sigma22 <- Sigma[(p1+1):p,(p1+1):p]
    x2 <- mu2 + Sigma21%*%solve(Sigma11,x1-mu1)
    return(x2)
}
M.estimate <- colMeans(train[,52:102])
SError <- CSCSError <- HuangError <- ShojError <- matrix(0,nrow=nrow(test),ncol=51)
for(k in 1:nrow(test)){
    if(trainsize>102){
        S.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=S)
        CSCS.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=CSCS)
        Huang.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=Huang)
        Shoj.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=Shoj)
        SError[k,] <- (S.mu2-as.numeric(test[k,52:102]))
        CSCSError[k,] <- (CSCS.mu2-as.numeric(test[k,52:102]))
        HuangError[k,] <- (Huang.mu2-as.numeric(test[k,52:102]))
        ShojError[k,] <- (Shoj.mu2-as.numeric(test[k,52:102]))
    } else {
        CSCS.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=CSCS)
        Huang.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=Huang)
        Shoj.mu2 <- predict.mean(test[k,1:51],mu=mu,Sigma=Shoj)
        CSCSError[k,] <- (CSCS.mu2-as.numeric(test[k,52:102]))
        HuangError[k,] <- (Huang.mu2-as.numeric(test[k,52:102]))
        ShojError[k,] <- (Shoj.mu2-as.numeric(test[k,52:102]))
    }
    
}

if(trainsize>102){
    E.S <- colMeans(abs(SError))
    E.CSCS <- colMeans(abs(CSCSError))
    E.Huang <- colMeans(abs(HuangError))
    E.Shoj <- colMeans(abs(ShojError))
} else {
    E.CSCS <- colMeans(abs(CSCSError))
    E.Huang <- colMeans(abs(HuangError))
    E.Shoj <- colMeans(abs(ShojError)) 
}

Time = 52:102
AE = E.S
Method = rep("S",length(AE))
Sdata = data.frame(Method,AE,Time)

AE = E.CSCS
Method = rep("CSCS",length(AE))
cscsdata = data.frame(Method,AE,Time)

AE = E.Huang
Method = rep("Sparse Cholesky",length(AE))
huangdata = data.frame(Method,AE,Time)

AE = E.Shoj
Method = rep("Sparse DAG",length(AE))
shojdata = data.frame(Method,AE,Time)

if(trainsize>102){
    mergeddata = rbind(Sdata,cscsdata,huangdata,shojdata)
    mergeddata = as.data.frame(mergeddata)
} else {
    mergeddata = rbind(cscsdata,huangdata,shojdata)
    mergeddata = as.data.frame(mergeddata)
}

pic<- ggplot(mergeddata, aes(x = Time, y=AE, group=Method, color=Method)) + 
    geom_line() 
### +theme(legend.position="none")
print(pic)

filename = paste(ord,trainsize,'callcenter-cv.pdf',sep="")
ggsave(filename,width = 8, height = 4)

#aafe
write(by(mergeddata$AE,(mergeddata$Method),sum), file = paste(trainsize,"aafe.txt",sep=""))

aafe <- by(mergeddata$AE,(mergeddata$Method),sum)

#min counts
write(table(aggregate(mergeddata$AE,list(mergeddata$Time),which.min)$x),file = paste(trainsize,"mincount.txt",sep=""))

min.count <- table(aggregate(mergeddata$AE,list(mergeddata$Time),which.min)$x)

if(trainsize>102){
    ll <- c(sll,cscsll,huangll,shojll)
    names(ll) <- c("S","CSCS","Sparse Cholesky", "Sparse Graph")
} else {
    ll <- c(cscsll,huangll,shojll)
    names(ll) <- c("CSCS","Sparse Cholesky", "Sparse Graph")
}

write(ll,file = paste(trainsize,"ll.txt",sep=""))
save(aafe,min.count,ll,file = paste(trainsize,".RData",sep=""))

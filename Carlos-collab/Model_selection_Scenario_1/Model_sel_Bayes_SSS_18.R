
library(Matrix)
library(BDgraph)
#setwd("WHATEVER")


#####################################################################################################
#####################################################################################################
########------------------------Program to read data from QMSim-----------------------------#########
#####################################################################################################
#####################################################################################################

#--To read genotypes
#4=1,2
QTL1=read.table("Pop1_qtl_004.txt",header=F)                                     
QTL1=as.matrix(QTL1[,-1])
QTL.Matrix=function(QTL){
	template=seq(1,ncol(QTL)-1,2)
	A=matrix(0,nrow=nrow(QTL),ncol=ncol(QTL))
	for(i in 1:nrow(QTL)){
		for(j in template){
			A[i,j]=QTL[i,j]+QTL[i,j+1]-3
			        }
			     }
	A=A[,template]
	return(A)
				}




###---Reading phenotypes and total additive effects     
#4=1,2                                   
data1=read.table("Pop1_data_004.txt",header=T)                                     
#var(cbind(data1$Res,data1$QTL,data1$Phen))
W=QTL.Matrix(QTL=QTL1)
nloci=ncol(W)


#########################################################
###-----Create training and testing population-----######
#########################################################

data1.Training=data1[-which(data1$G==3), ]
W.Train=W[-which(data1$G==3), ]
n.train=nrow(W.Train)
data1.Test=data1[which(data1$G==3),]


####---Reading graph files----#####

Graph4_9=read.table("Graph4 9 .txt",header=F)


###############################################################
########---Functions to build the adjacency matrices----#######
###############################################################

Adj=function(Nodes,nloci){
	library(Matrix)
	Adj=matrix(0,nrow=nloci,ncol=nloci)
	for(i in 1:nrow(Nodes)){
		Adj[Nodes[i,1],Nodes[i,2]]=1
                              }
	return(Adj)
			       }


#---This second version builds a symmetric matrix instead of an upper triangular one

Adj2=function(Nodes,nloci){
	library(Matrix)
	Adj=diag(nloci)
      for(i in 1:nrow(Nodes)){
		Adj[Nodes[i,1],Nodes[i,2]]=1
            Adj[Nodes[i,2],Nodes[i,1]]=1
                              }
	return(Adj)
			       }


#########----Function to do full MC integration (i.e.,NO Laplace approximation)-------#########

MCint2=function(y,W,Tau=210,V=4.1,nsamples=2000,Adjac){
 nloci=ncol(W)
  library(mvtnorm)
   library(pscl)
    n=length(y)
     I=diag(n)
      WTW=crossprod(W)
	 yt=t(y)
	  Wt=t(W)
	   funct=matrix(0,nrow=nsamples,ncol=1)
	    for(i in 1:nsamples){
         resvar=rigamma(1,alpha=V/2,beta=Tau/2)
        Omega=matrix(round(rgwish(n=1,adj.g=Adjac,b=10,D=diag(nloci)),8),ncol=nloci,byrow=TRUE)
       Sigma=as.matrix(forceSymmetric(chol2inv(chol(Omega))))
      g=matrix(rmvnorm(1,mean=matrix(0,nrow=nloci,ncol=1),sigma=Sigma,method="svd"))
     #funct[i]=dmvnorm(yt,mean=W%*%g,sigma=diag(resvar,nrow=n),log=TRUE)
     funct[i]=dmvnorm(yt,mean=W%*%g,sigma=diag(resvar,nrow=n),log=FALSE)
    }
   estimate=mean(funct)
  var.estimate=var(funct)*(nsamples-1)/(nsamples^2)
 return(list(estimate,var.estimate,funct))        
}




#####################################################################################################################
#####################################################################################################################
######-----Syed's function to perform the SSS algorithm of Ben-David et al. (2015)--------###########################
#####################################################################################################################
#####################################################################################################################


getNgraphs <- function(D,N){
    N1 <- list()
    dimD <- dim(D)
    for(i in 1:N){
        N1[[i]] <- D
        x <- floor(runif(1,1,dimD[1]))
        y <- floor(runif(1,1,x))
        N1[[i]][x,y] = ifelse(N1[[i]][x,y]==0,1,0)
    }
    N1
}

 
 
getscores <- function(graphs){
    scores <- rep(0,length(graphs))
    for(i in 1:length(graphs)){
        scores[i] <-100*MCint2(y=data1.Training$Phen,W=W.Train,Tau=210,V=4.1,nsamples=2000,Adjac=graphs[[i]])[[1]]
        #scores[i] <- gnorm(adj.g=graphs[[i]],b=10,D=diag(nloci),iter=1000) ###This one is for m>n case only
        }
    scores
}

 

 
getDnew <- function(gamma, scores2, graphs2){
    pvec <- exp(scores2^gamma)
    pvec <- pvec/sum(pvec)
    cpvec <- cumsum(pvec)
    u <- runif(1,0,1)
    Dind <- min(sum(u > cpvec) + 1,length(scores2))
    graphs2[[Dind]]
}

 
getLk <- function(gamma,M,N1,D0){
 
    graphs <- getNgraphs(D0,N1)
    scores <- getscores(graphs)
    D0 <- getDnew(gamma,scores,graphs)
    for(i in 2:M){
        graphs2 <- getNgraphs(D0,N1)
        scores2 <- getscores(graphs2)
        D0 <- getDnew(gamma,scores2,graphs2)
        graphs <- append(graphs,graphs2)
        scores <- append(scores,scores2)
    }
    return(list(graphs=graphs, scores=scores))
}

 
 
Largestscoregraph <- function(gamma,D,M,N1){
    L <- list()
    for(i in 1:length(D)){
        newL <- getLk(gamma,M,N1,D[[i]])
        L$graphs <- append(newL$graphs,L$graphs)
        L$scores <- append(newL$scores,L$scores)
    }
    L$graphs[[which.max(L$scores)]]
}

####---Example to perform the SSS algorithm for a particular k----####

Adjac=Adj(Nodes=Graph4_9,nloci=nloci)

Score4_9=getLk(gamma=0.5,M=3,N1=10,D0=Adjac)

save(Score4_9,file=paste("score",toString(4),toString(9),".rdata",sep=""))

####--So, this gives L(k) for a given k, k=1,2,...,15. Instead of using Largestscoregraph, 
####--we run different batches of graphs in parallel, this will be faster.




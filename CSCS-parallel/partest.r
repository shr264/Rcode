library(Matrix)
library(MASS)
library(microbenchmark)
source("CSCS2.r")
source("parCSCS.r")


p = 100  # predictors
n = 500  ############# this needs to change to 250
z = 0.3  # fraction of non-zeros in cholesky paramter of L
s = 0.5  # fraction of negative coefficients
a = 0.3  # minimum magnitude of non-zero coefficients
b = 0.7  # maximum magnitude of non-zero coefficients

plower = p*(p-1)/2

set.seed(12345) #seed for generating L
## diagonals
D = runif(p,2,5)

## off-diagonals
T = diag(p)
T[upper.tri(T)] = 0
T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) * 
		   ifelse(runif(plower)<z,  1, 0) * 
		   runif(plower, a, b))

L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
omega = t(L) %*% L            # omega
sigma = solve(omega)
sigma[abs(sigma)<1e-10] = 0   # set numerical error to zero


sum(abs(omega)>0)/choose(p,2) #sparsity in omega

set.seed(23456 + 1) #seed for generating data
X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations

X = scale(X, center = TRUE, scale = TRUE) # centered obs

microbenchmark(CSCS2(Y=X,lambda=0.5),parCSCS(Y=X,lambda=0.5,clusters = 3), times=5)


registerDoParallel(cores=2)
out = foreach(i=1:4) %dopar% {a = 1:i
                        a}

as.matrix(unlist(out), byrow = TRUE)

L = diag(4)
L[lower.tri(L,diag=TRUE), byrow] = unlist(out)



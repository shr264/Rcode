CSCS2 <- function( Y, lambda, L=NULL, maxitr=100, tol=1e-4, warmstart=FALSE) {

  #### inputs
  ## Y: n by p matrix of data
  ## lambda: l1-penalty parameter
  ## maxitr: maximum number of iterations allowed (diagonal/off-diagonal optimization)
  ## tol: if maximum difference of iterates is less than tol, consider converged
  ## warmstart: warmstarting actually made the runs slower (overhead may be too expensive)

  #### outputs
  ## L: lower triangular matrix of cholesky factor computed with CSCS algorithm
  ## itr_log: (p-1)-vector of number of iterations 
  ## eps_log: (p-1)-vector of number maximum difference for considering convergence
  
  
  n = nrow(Y)
  p = ncol(Y)
  
  if (length(lambda)==1) lambda = rep(lambda,p)

  if (is.null(L)) L = diag(p)

  S = (t(Y)%*%Y)/n

  itr_log = eps_log = NULL

  L[1, 1] = 1/sqrt(S[1,1])
  for (k in 2:p){ ## Loop in Algorithm 2

    ## nu_k vector (equation 2.5)
    nuk_old = nuk_new = c(rep(0, k-1), 1) 
    r = 0

    repeat {      ## Loop in Algorithm 1

      r = r + 1
      
      km1_ = 1:(k-1)    ## 1, ..., k-1 is off-diagonal elements indices

      ## Update off-diagonal terms
      hk = lassoshooting(XtX    =  S[km1_, km1_, drop=FALSE], 
			 Xty    = -nuk_old[k] * S[km1_, k],
			 lambda =  0.5*lambda[k])
      nuk_new[km1_] = hk$coefficients

      ## Update diagonal term
      sumterm = sum(nuk_new[km1_] * S[k, km1_])
      nuk_new[k] = (-sumterm + sqrt(sumterm^2 + 4*S[k, k]))/(2*S[k, k])

      ## Check convergence
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
	L[k, 1:k] = nuk_new
	eps_log = c(eps_log, maxdiff)
	itr_log = c(itr_log, r)
	break
      } else {
	nuk_old = nuk_new
      }
    }
  }

  list(L=L, itr=itr_log, eps=eps_log)

}

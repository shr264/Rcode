huang <- function( Y, lambda, maxitr=100, tol=1e-4, returnL=FALSE){

  #### inputs
  ## Y: n by p matrix of data
  ## lambda: l1-penalty parameter
  ## maxitr: maximum number of iterations allowed (diagonal/off-diagonal optimization)
  ## tol: if maximum difference of iterates is less than tol, consider converged
  ## returnL: huang algorithm returns T and D; but, instead, compute and return L 

  #### outputs
  ## T: lower triangular matrix of correlation coefficients
  ## D: p-vector of marginal variances
  ## L: lower triangular matrix of cholesky factor computed with CSCS algorithm
  ## itr_log: (p-1)-vector of number of iterations 
  ## eps_log: (p-1)-vector of number maximum difference for considering convergence
  
  p = ncol(Y)
  n = nrow(Y)
  
  if (length(lambda)==1) lambda = rep(lambda,p)

  S = (t(Y) %*% Y)/n

  T = diag(p)      # identity matrix of size p
  D = rep(0, p)    # vector of zeros
  D[1] = D[1] = S[1,1] # D_11 = S_11

  itr_log = eps_log = NULL

  for (k in 2:p){

    ## nu_k[1:(k-1)] = phi_ij and nu_k[k] = sigma_k (or D_k)
    ##  (equations 4 and 5)
    nuk_old = nuk_new = c(rep(0, k-1), 1)
    r = 0
    km1_ = 1:(k-1)

    repeat {

      r = r + 1

      ## Update diagonals
      nuk_new[k] = S[k,k] + 2*sum(S[k,km1_] * nuk_old[km1_]) + 
	quad.form(S[km1_,km1_], nuk_old[km1_])

      ## Update off-diagonals
      output = lassoshooting(XtX= S[km1_,km1_,drop=FALSE], 
			      Xty=-S[km1_,k], lambda=0.5*nuk_new[k]*lambda[k])
      nuk_new[km1_] = output$coefficients

      ## Check convergence
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
	T[k, km1_] = nuk_new[km1_]
	D[k] = nuk_new[k]
	eps_log = c(eps_log, maxdiff)
	itr_log = c(itr_log, r)
	break
      } else {
	nuk_old = nuk_new
      }

    }

  }

  if (returnL){
    output = list(L=T%*%diag(1/sqrt(D)), itr=itr_log, eps=eps_log)
  } else {
    output = list(T=T, D=D,              itr=itr_log, eps=eps_log)
  }

  output

}

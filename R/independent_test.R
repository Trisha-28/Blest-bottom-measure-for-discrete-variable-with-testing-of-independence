#'@title  Testing independence based on Blest measure  between random variables
#'
#'@description This function computes the matrix of Blest statistics and p-value between random variables with discrete distributions.
#'
#'@param x           Data matrix
#'
#'@return \item{stat}{Blest  statistics from multilinear copula}
#'@return \item{pvalue}{P-values for the test statistics}
#'
#'x <-rpois(100,3)
#'y <-rpois(100,4)
#'z <-cbind(x,y)
#'out<-EstDep(z)
#'
#'
EstDep = function(x)
{
  dd = dim(x)
  n = dd[1]
  d = 2
  blest = matrix(0,ncol=2,nrow=2)
  pvalblest = matrix(0,ncol=2,nrow=2)
  blest[1,1] = 1
  blest[2,2] = 1
  out0 = stat_dep(x[,1],x[,2])
  out1 = stat_dep(x[,2],x[,1])
  blest[1,2] = out0$blest
  blest[2,1] = out1$blest
  zblest_ij = sqrt(n) * blest[1,2]
  zblest_ji = sqrt(n) * blest[2,1]
  pvalblest[1,2] = 200 * pnorm(-abs(zblest_ij))
  pvalblest[2,1] = 200 * pnorm(-abs(zblest_ji))
  pvalblest[1,1] = 100
  pvalblest[2,2] = 100
  out = list(blest = blest, pvalblest = pvalblest)
  return(out)
}

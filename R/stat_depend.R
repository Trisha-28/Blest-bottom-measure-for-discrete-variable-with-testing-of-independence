#'@title Computes the Blest measure for tests of independence
#'
#'@description This function Computes the Blest measure for tests of independence used in indep_test
#'@param x     Vector of length n
#'@param y     Vector of length n
#'
#'@return \item{rho}{Blest measure}
#'
#'@examples
#' x <-rpois(100,3)
#' y <-rpois(100,4)
#' out <- stat_dep(x,y)

#'
#'@keywords internal
#'
#'@export
#'
#'


stat_dep = function(x, y){
  n = length(x)

  out = .C("estdep",
           as.double(x),
           as.double(y),
           as.integer(n),
           blest = double(1),
           std = double(1),
           PACKAGE = "Blestbottom"
  )

  return(list(
    blest = out$blest
  ))
}


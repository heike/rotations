arsample <- function(f,g,M, kappa, ...) {
  found=FALSE
  while(!found) {
    x <- g(1, ...)
    y <- runif(1, min=0, max=M)
    if (y < f(x, kappa)) found=TRUE
  }
  return(x)
}

#' Random simulation from angle distributions 
#'
#' acceptance-rejection random sampling from angle distributions
#'
#' @export
#' @param n sample size
#' @param f density to be sampled from
#' @param g density, enveloping f, i.e. f < g 
#' @param M real valued constant
#' @param ... parameters for densities f and g
#' @examples
#' # sample from haar distribution
#' x <- rar(10000, haar, runif, 1/pi, min=-pi, max=pi)
#' 
#' kappa=0.5
#' M <- max(fisher(seq(-pi, pi, length=1000), kappa))
#' x.fisher <- rar(10000, fisher, runif, M, min=-pi, max=pi, kappa=kappa)
rar <- function(n, f,g, M, ...) {
  res <- vector("numeric", length=n)
  for (i in 1:n) res[i] <- arsample(f,g,M, ...)
  return(res)
}



#' Angle distributions 
#'
#' Densities of the most commonly used angle distributions: Haar, Fisher, Cayley, and von Mises Distributions
#'
#' @export
#' @param x vector of angles (between -pi and pi)
#' @param kappa concentration parameter
haar <- function(x) return((1-cos(x))/(2*pi))
haarkappa <- function(x, kappa=0) return(haar(x))


fisher <- function(x, kappa) {
	exp(2*kappa*cos(x))*(1-cos(x))/(besselI(kappa,0)-besselI(kappa,1))/(2*pi)
}
cayley <- function(x, kappa) {
	return(.5/sqrt(pi)*gamma(kappa+2)/gamma(kappa+0.5)*2^(-kappa)*(1+cos(x))^kappa*(1-cos(x)))
}
mises <- function(x, kappa) return(exp(kappa*cos(x))/(2*pi*besselI(kappa,0)))

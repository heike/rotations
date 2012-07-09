#' Angle distributions 
#'
#' Densities of the most commonly used angle distributions: Haar, Fisher, Cayley, and von Mises.
#'
#' The Haar density for angles is given as \eqn{(1-\cos(x))/(2*\pi)}, where \eqn{x \in (-\pi, \pi)}. When the fix axis \eqn{U} of rotation \eqn{R} is randomly picked, rotations with angles from the Haar measure are uniformly distributed on the space of all 3d rotations, \eqn{SO(3)}.
#' @export
#' @param x vector of angles (between -pi and pi)
#' @param kappa concentration parameter, kappa > 0
#' @aliases fisher cayley mises
#' @usage haar <- function(x) 
#'   fisher <- function(x, kappa) 
#'   cayley <- function(x, kappa)
#'   mises <- function(x, kappa)
#' @examples
#' x <- seq(-pi, pi, by=0.01)
#' plot(x, haar(x), type="l")
#' lines(x, fisher(x, kappa=1), col=2)
#' lines(x, cayley(x, kappa=1), col=3)
#' 
#' plot(x, fisher(x, kappa=1)/haar(x), type="l")
#' lines(x, cayley(x, kappa=1)/haar(x), col=3)
#' # Cayley with the same concentration parameter is heavier tailed than Fisher. 
haar <- function(x) return((1-cos(x))/(2*pi))
haarkappa <- function(x, kappa=0) return(haar(x))

#' @export
fisher <- function(x, kappa) {
	exp(2*kappa*cos(x))*(1-cos(x))/(besselI(kappa,0)-besselI(kappa,1))/(2*pi)
}

#' @export
cayley <- function(x, kappa) {
	return(.5/sqrt(pi)*gamma(kappa+2)/gamma(kappa+0.5)*2^(-kappa)*(1+cos(x))^kappa*(1-cos(x)))
}

#' @export
mises <- function(x, kappa) return(exp(kappa*cos(x))/(2*pi*besselI(kappa,0)))


#' Angular distributions in the rotations packages
#' 
#' Density and random variate generation for symmetric probability distributions in the rotations package
#' 
#' The functions for the density function and random variate generation are named in the usual form dxxxx and rxxxx 
#' respectively.        
#' \itemize{
#' 	\item See \code{\link{dcayley}} for the Cayley distribution.
#' 	\item See \code{\link{dfisher}}for the matrix Fisher distribution.
#' 	\item See \code{\link{dhaar}} for the uniform distribution on the circle.
#' 	\item See \code{\link{dvmises}} for the von Mises-Fisher distribution.
#' }
#' 
#' @name Angular-distributions

NULL

#' Sample of size n from target density f
#'
#' @author Heike Hofmann
#' @param n number of sample wanted
#' @param f target density
#' @param M maximum number in uniform proposal density
#' @param ... additional arguments sent to arsample
#' @return a vector of size n of observations from target density

rar <- function(n, f, M, ...) {
  res <- vector("numeric", length = n)
  for (i in 1:n) res[i] <- arsample.unif(f, M, ...)
  return(res)
}


#' The Symmetric Cayley Distribution
#'
#' Density and random generation for the Cayley distribution with concentration kappa
#'
#' The symmetric Cayley distribution with concentration kappa (or circular variance nu) had density 
#' \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r).}
#'
#' @name Cayley
#' @aliases Cayley rcayley dcayley
#' @usage dcayley(r, kappa = 1, nu = NULL, Haar = TRUE, lower.tail=TRUE)
#' @usage rcayley(n, kappa = 1, nu = NULL)
#' @param r vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa Concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return \code{dcayley} gives the density, \code{rcayley} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package
#' @cite Schaeben97 leon06

NULL


#' @rdname Cayley
#' @aliases Cayley rcayley dcayley
#' @export

dcayley <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
 	den <- 0.5 * gamma(kappa + 2)/(sqrt(pi) * 2^kappa * gamma(kappa + 0.5)) * (1 + cos(r))^kappa * (1 - cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
}


#' @rdname Cayley
#' @aliases Cayley rcayley dcayley
#' @export

rcayley <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  bet <- rbeta(n, kappa + 0.5, 3/2)
  theta <- acos(2 * bet - 1) * (1 - 2 * rbinom(n, 1, 0.5))
  return(theta)
}

#' The Matrix Fisher Distribution
#'
#' Density and random generation for the matrix Fisher distribution with concentration kappa
#'
#' The matrix Fisher distribution with concentration kappa (or circular variance nu) has density
#' \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  
#' \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#'
#' @name Fisher
#' @aliases Fisher dfisher rfisher
#' @usage dfisher(r, kappa = 1, nu = NULL, Haar = TRUE, lower.tail=TRUE)
#' @usage rfisher(n, kappa = 1, nu = NULL)
#' @param r vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa concentration paramter
#' @param nu circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return \code{dfisher} gives the density, \code{rfisher} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Fisher
#' @aliases Fisher dfisher rfisher
#' @export

dfisher <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  n<-length(r)
  den<-rep(0,n)
  
 	den <- exp(2 * kappa * cos(r)) * (1 - cos(r))/(2 * pi * (besselI(2 * kappa, 0) - besselI(2 * kappa, 1)))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
  
}

#' @rdname Fisher
#' @aliases Fisher dfisher rfisher
#' @export


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  M <- max(dfisher(seq(-pi, pi, length = 1000), kappa,Haar=F))
  return(rar(n, dfisher, M, kappa = kappa, Haar=F))
}

#' Haar Measure
#'
#' Uniform density on the circle
#' 
#' The uniform density on the circle  (also referred to as Haar measure)
#' has the density \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}.}
#'
#' @name Haar
#' @aliases Haar dhaar rhaar
#' @usage dhaar(r, lower.tail=TRUE)
#' @usage rhaar(n)
#' @param r Where the density is being evaluated
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @return \code{dhaar} gives the density, \code{rhaar} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Haar
#' @aliases Haar dhaar rhaar
#' @export

dhaar <- function(r){
	
	den <-(1 - cos(r))/(2 * pi)
	
	return(den)
} 

#' @rdname Haar
#' @aliases Haar dhaar rhaar
#' @export


rhaar<-function(n){
	
	lenn<-length(n)
	if(lenn>1)
		n<-lenn
	
  return(rar(n, dhaar, 1/pi))
}

#' The circular-von Mises distribution
#'
#' Density for the the circular von Mises-based distribution with concentration kappa
#' 
#' The circular von Mises-based distribution has the density
#' \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa cos(r)}.}
#'
#' @name Mises
#' @aliases Mises dvmises rvmises
#' @usage dvmises(r, kappa = 1, nu = NULL, Haar = TRUE, lower.tail=TRUE)
#' @usage rvmises(n, kappa = 1, nu = NULL)
#' @param r vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return \code{dvmises} gives the density, \code{rvmises} generates random deviates
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Mises
#' @aliases Mises dvmises rvmises
#' @export

dvmises <- function(r, kappa = 1, nu = NULL, Haar = T) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  den <- 1/(2 * pi * besselI(kappa, 0)) * exp(kappa * cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
}

#' @rdname Mises
#' @aliases Mises dvmises rvmises
#' @export

rvmises <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  u <- runif(3, 0, 1)
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a))/(2 * kappa)
  r <- (1 + b^2)/(2 * b)
  theta <- rep(10, n)
  
  for (i in 1:n) {
    
    while (theta[i] == 10) {
      # Step 1
      u <- runif(3, 0, 1)
      z <- cos(pi * u[1])
      f <- (1 + r * z)/(r + z)
      c <- kappa * (r - f)
      
      # Step 2
      u <- runif(3, 0, 1)
      if ((c * (2 - c) - u[2]) > 0) {
        
        theta[i] = sign(u[3] - 0.5) * acos(f)
        
      } else {
        
        if (log(c/u[2]) + 1 - c < 0) {
          u <- runif(3, 0, 1)
        } else {
          u <- runif(3, 0, 1)
          theta[i] = sign(u[3] - 0.5) * acos(f)
        }
      }
    }
  }
  return(theta)
}



#' UARS density function
#' 
#' Evaluate the UARS density with a given angular distribution.
#' 
#' @param o Value at which to evaluate the UARS density
#' @param S principal direction of the distribution
#' @param kappa concentration of the distribution
#' @param c Circular distribution, choices: 'dcayley', 'dfisher', 'dhaar' and 'dvmises'
#' @return density value at o
#' @export

duars<-function(o,S=diag(3),kappa,c,...){
	
	o<-matrix(o,3,3)
	trStO<-sum(diag(t(S)%*%o))
	r<-acos(.5*(trStO-1))
	cr<-c(r,kappa,...)
	den<-4*pi*cr/(3-trStO)
	return(den)
}
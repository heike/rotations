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
#' Density of the symmetric Cayley distribution with concentration kappa
#'
#' The symmetric Cayley distribution with concentration kappa (or circular variance nu) had density 
#' \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r).}
#'
#' @param r vector of quantiles
#' @param kappa Concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @param lower.tail logical; if TRUE probabilites are \eqn{P(X\le x)}
#' @return Cayley density with concentration kappa evaluated at r
#' @seealso \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dvmises}}
#' @export
#' @cite Schaeben97 leon06

dcayley <- function(r, kappa = 1, nu = NULL, Haar = TRUE, lower.tail=TRUE) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
 	den <- 0.5 * gamma(kappa + 2)/(sqrt(pi) * 2^kappa * gamma(kappa + 0.5)) * (1 + cos(r))^kappa * (1 - cos(r))
  
  if(!lower.tail)
  	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
}

#' The Symmetric Cayley Distribution
#'
#' Random generation from the symmetric Cayley distribution with concentration kappa (or circular variance nu)
#' 
#' The symmetric Cayley distribution with concentration kappa (or circular variance nu) has density 
#' \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r).}
#'
#' @param n sample size
#' @param kappa The concentration paramter
#' @param nu An alternative to kappa; circular variance
#' @return random deviates from Cayley distribution with concentration kappa
#' @cite Schaeben97 leon06
#' @seealso \code{\link{dcayley}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}
#' @export
#' @examples
#' r<-rcayley(20,0.01)

rcayley <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
  bet <- rbeta(n, kappa + 0.5, 3/2)
  theta <- acos(2 * bet - 1) * (1 - 2 * rbinom(n, 1, 0.5))
  return(theta)
}

#' The Matrix Fisher Distribution
#'
#' Density of the matrix Fisher distribution with concentration kappa
#'
#' The matrix Fisher distribution with concentration kappa (or circular variance nu) has density
#' \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  
#' \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#'
#' @param r vector of quantiles
#' @param kappa concentration paramter
#' @param nu circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @param lower.tail logical; if TRUE probabilites are \eqn{P(X\le x)}
#' @return value of Fisher matrix distribution with concentration \eqn{\kappa} evaluated at r
#' @seealso \code{\link{rfisher}}, \code{\link{dhaar}},\code{\link{dvmises}},\code{\link{dcayley}}
#' @export

dfisher <- function(r, kappa = 1, nu = NULL, Haar = TRUE, lower.tail=TRUE) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  n<-length(r)
  den<-rep(0,n)
  
  for(i in 1:n)
  	den[i] <- exp(2 * kappa * cos(r[i])) * (1 - cos(r))/(2 * pi * (besselI(2 * kappa, 0) - besselI(2 * kappa, 1)))
  
  if(!lower.tail)
  	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
  
}

#' The Matrix Fisher Distribution
#'
#' Random generation for the matrix Fisher distribution with concentration kappa (or circular variance nu)
#' 
#' The matrix Fisher distribution with concentration kappa (or circular variance nu) has density
#' \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  
#' \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#'
#' @param n sample size
#' @param kappa the concentration parameter
#' @param nu An alternative to kappa; circular variance
#' @return a sample of size \eqn{n} from the matrix Fisher distribution with concentration \eqn{\kappa}
#' @seealso \code{\link{dfisher}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}
#' @export


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
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
#' @param r Where the density is being evaluated
#' @param lower.tail logical; if TRUE probabilites are \eqn{P(X\le x)}
#' @return the probability density evaluated at r
#' @seealso \code{\link{rhaar}}, \code{\link{dfisher}},\code{\link{dvmises}},\code{\link{dcayley}}
#' @export

dhaar <- function(r, lower.tail = TRUE){
	
	den <-(1 - cos(r))/(2 * pi)
	
	if(lower.tail) return(den) else return(1-den)
} 

#' Haar Measure
#'
#' Random generation from the circle
#' 
#' The uniform density on the circle  (also referred to as Haar measure)
#' has the density \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}.}
#' 
#' @param n sample size
#' @return a sample of size \eqn{n} from the uniform distribution on the sphere
#' @seealso \code{\link{dhaar}},\code{\link{rfisher}},\code{\link{rvmises}},\code{\link{rcayley}}
#' @export


rhaar<-function(n){
  return(rar(n, dhaar, 1/pi))
}

#' The circular-von Mises distribution
#'
#' Density for the the circular von Mises-based distribution with concentration kappa
#' 
#' The circular von Mises-based distribution has the density
#' \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.
#'
#' @param r vector of quantiles
#' @param kappa concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @param lower.tail logical; if TRUE probabilites are \eqn{P(X\le x)}
#' @export
#' @return value of circular-von Mises distribution with concentration \eqn{\kappa} evaluated at r
#' @seealso \code{\link{rvmises}}, \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dcayley}}

dvmises <- function(r, kappa = 1, nu = NULL, Haar = T, lower.tail=TRUE) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  den <- 1/(2 * pi * besselI(kappa, 0)) * exp(kappa * cos(r))
  
  if(!lower.tail)
  	den<-1-den
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
}

#' The circular-von Mises distribution
#'
#' Random generation for the the circular von Mises-based distribution with concentration kappa
#' 
#' The circular von Mises-based distribution has the density
#' \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.
#'
#' @param n The number of angles desired
#' @param kappa The concentration parameter of the distribution
#' @param nu An alternative to kappa; circular variance
#' @return S3 \code{rvmises} object; a vector of n angles following the von Mises Circular distribution with concentration kappa and mean/mode 0
#' @export
#' @examples
#' r<-rvmises(20,0.01)

rvmises <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
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

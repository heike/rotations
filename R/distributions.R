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


#' Simulate misorientation angles from Cayley distribution
#'
#' This function allows the user to simulate \eqn{n} misorientation angles from the Cayley distribution symmetric about 0 on interval \eqn{(-\pi,\pi]}.  The relationship between Cayley and Beta distribution is used.
#' The symmetric Cayley distribution has a density of the form \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r)}.
#' It was orignally given in the material sciences literature by Schaben 1997 and called the de la Vallee Poussin distribution but was more recently discussed and
#' introduced in a more general manner by Leon 06.
#'
#' @param n sample size
#' @param kappa The concentration paramter
#' @param nu An alternative to kappa; circular variance
#' @return vector of n observations from Cayley(kappa) distribution
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

#' Simulate a data set of size \eqn{n} from the matrix Fisher angular distribution
#'
#' The symmetric matrix fisher distribution has the density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#' This function allows for simulation of \eqn{n} random deviates with density \eqn{C_\mathrm{{F}}(r|\kappa)} and \eqn{\kappa} provided by the user.
#'
#' @param n sample size
#' @param kappa the concentration parameter
#' @param nu An alternative to kappa; circular variance
#' @return a sample of size \eqn{n} from the matrix Fisher distribution with concentration \eqn{\kappa}
#' @seealso \code{\link{dfisher}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  M <- max(dfisher(seq(-pi, pi, length = 1000), kappa,Haar=F))
  return(rar(n, dfisher, M, kappa = kappa, Haar=F))
}

#' Simulate a data set of size \eqn{n} from the uniform distribution on the sphere
#'
#' The uniform distribution has the density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi}1-\cos(r)}.  The is also 
#' know as the Haar measure.
#' 
#' @param n sample size
#' @return a sample of size \eqn{n} from the uniform distribution on the sphere
#' @seealso \code{\link{dhaar}},\code{\link{rfisher}},\code{\link{rvmises}},\code{\link{rcayley}}


rhaar<-function(n){
  return(rar(n, dhaar, 1/pi))
}


#' Generate a vector of angles(r) from the von Mises Circular distribution
#'
#' The circular von Mises-based distribution has the density \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.  This function allows the use to
#' simulate \eqn{n} random deviates from \eqn{C_\mathrm{M}(r|\kappa)} given a concentration parameter \eqn{\kappa}.
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

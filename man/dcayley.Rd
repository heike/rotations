\name{dcayley}
\alias{dcayley}
\title{Symmetric Cayley distribution for angular data}
\usage{
  dcayley(r, kappa = 1, nu = NULL, Haar = T)
}
\arguments{
  \item{r}{Where the density is being evaluated}

  \item{kappa}{The concentration paramter, taken to be
  zero}

  \item{nu}{The circular variance, can be used in place of
  kappa}

  \item{Haar}{logical, if density is evaluated with respect
  to Haar measure or Lebesgue}
}
\value{
  value of Cayley distribution with concentration
  \eqn{\kappa} evaluated at r
}
\description{
  The symmetric Cayley distribution has a density of the
  form \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}}
  \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos
  r)^\kappa(1-\cos r)}. It was orignally given in the
  material sciences literature by Schaben 1997 and called
  the de la Vallee Poussin distribution but was more
  recently discussed and introduced in a more general
  manner by Leon 06.
}
\references{
  Schaeben H (1997). "A Simple Standard Orientation Density
  Function: The Hyperspherical de la Vallée Poussin
  Kernel." _Phys. Stat. Sol. (B)_, *200*, pp. 367-376.

  Le\'on C, é JM and Rivest L (2006). "A statistical model
  for random rotations." _Journal of Multivariate
  Analysis_, *97*(2), pp. 412-430.
}
\seealso{
  \code{\link{rcayley}},\code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dvmises}}
}


\name{dcayley}
\alias{dcayley}
\title{The Symmetric Cayley Distribution}
\usage{
  dcayley(r, kappa = 1, nu = NULL, Haar = TRUE,
    lower.tail = TRUE)
}
\arguments{
  \item{r}{vector of quantiles}

  \item{kappa}{Concentration paramter}

  \item{nu}{The circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}

  \item{lower.tail}{logical; if TRUE probabilites are
  \eqn{P(X\le x)}}
}
\value{
  Cayley density with concentration kappa evaluated at r
}
\description{
  Density of the symmetric Cayley distribution with
  concentration kappa
}
\details{
  The symmetric Cayley distribution with concentration
  kappa (or circular variance nu) had density
  \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}}
  \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos
  r)^\kappa(1-\cos r).}
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
  \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dvmises}}
}


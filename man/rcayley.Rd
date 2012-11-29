\name{rcayley}
\alias{rcayley}
\title{The Symmetric Cayley Distribution}
\usage{
  rcayley(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{n}{sample size}

  \item{kappa}{The concentration paramter}

  \item{nu}{An alternative to kappa; circular variance}
}
\value{
  vector of n observations from Cayley(kappa) distribution
}
\description{
  This function allows the user to simulate \eqn{n}
  misorientation angles from the Cayley distribution
  symmetric about 0 on interval \eqn{(-\pi,\pi]}.  The
  relationship between Cayley and Beta distribution is
  used. The symmetric Cayley distribution has a density of
  the form \deqn{C_\mathrm{C}(r
  |\kappa)=\frac{1}{\sqrt{\pi}}
  \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos
  r)^\kappa(1-\cos r)}. It was orignally given in the
  material sciences literature by Schaben 1997 and called
  the de la Vallee Poussin distribution but was more
  recently discussed and introduced in a more general
  manner by Leon 06.
}
\examples{
r<-rcayley(20,0.01)
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
  \code{\link{dcayley}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}
}


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
  random deviates from Cayley distribution with
  concentration kappa
}
\description{
  Random generation from the symmetric Cayley distribution
  with concentration kappa (or circular variance nu)
}
\details{
  The symmetric Cayley distribution with concentration
  kappa (or circular variance nu) has density
  \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}}
  \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos
  r)^\kappa(1-\cos r).}
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


\name{rfisher}
\alias{rfisher}
\title{The Matrix Fisher Distribution}
\usage{
  rfisher(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{n}{sample size}

  \item{kappa}{the concentration parameter}

  \item{nu}{An alternative to kappa; circular variance}
}
\value{
  a sample of size \eqn{n} from the matrix Fisher
  distribution with concentration \eqn{\kappa}
}
\description{
  Random generation for the matrix Fisher distribution with
  concentration kappa (or circular variance nu)
}
\details{
  The matrix Fisher distribution with concentration kappa
  (or circular variance nu) has density
  \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
  where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel
  function of order \eqn{p} defined as
  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos
  r}dr}.
}
\seealso{
  \code{\link{dfisher}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}
}


\name{dfisher}
\alias{dfisher}
\title{The Matrix Fisher Distribution}
\usage{
  dfisher(r, kappa = 1, nu = NULL, Haar = TRUE,
    lower.tail = TRUE)
}
\arguments{
  \item{r}{vector of quantiles}

  \item{kappa}{concentration paramter}

  \item{nu}{circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}

  \item{lower.tail}{logical; if TRUE probabilites are
  \eqn{P(X\leq x)}}
}
\value{
  value of Fisher matrix distribution with concentration
  \eqn{\kappa} evaluated at r
}
\description{
  Density of the matrix Fisher distribution with
  concentration kappa
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
  \code{\link{rfisher}},
  \code{\link{dhaar}},\code{\link{dvmises}},\code{\link{dcayley}}
}


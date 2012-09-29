\name{dfisher}
\alias{dfisher}
\title{von Mises-Fisher distribution for angular data}
\usage{
  dfisher(r, kappa = 1, nu = NULL, Haar = T)
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
  value of Fisher matrix distribution with concentration
  \eqn{\kappa} evaluated at r
}
\description{
  The symmetric matrix fisher distribution has the
  density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
  where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel
  function of order \eqn{p} defined as
  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos
  r}dr}. This function allows the user to evaluate the
  function \eqn{C_\mathrm{{F}}(r|\kappa)} at \eqn{r} with
  \eqn{\kappa} provided by the user.
}
\seealso{
  \code{\link{rfisher}},
  \code{\link{dhaar}},\code{\link{dvmises}},\code{\link{dcayley}}
}


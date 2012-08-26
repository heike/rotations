\name{rhaar}
\alias{rhaar}
\title{Simulate a data set of size \eqn{n} from the uniform distribution on the sphere}
\usage{
  rhaar(n)
}
\arguments{
  \item{n}{sample size}
}
\value{
  a sample of size \eqn{n} from the uniform distribution on
  the sphere
}
\description{
  The uniform distribution has the
  density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi}1-\cos(r)}.
  The is also know as the Haar measure.
}
\seealso{
  \code{\link{dhaar}},\code{\link{rfisher}},\code{\link{rvmises}},\code{\link{rcayley}}
}


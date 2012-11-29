\name{rhaar}
\alias{rhaar}
\title{Haar Measure}
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
  Random generation from the circle
}
\details{
  The uniform density on the circle (also referred to as
  Haar measure) has the density
  \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}.}
}
\seealso{
  \code{\link{dhaar}},\code{\link{rfisher}},\code{\link{rvmises}},\code{\link{rcayley}}
}


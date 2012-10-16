\name{angle}
\alias{angle}
\alias{angle.EA}
\alias{angle.Q4}
\alias{angle.SO3}
\title{Find the angle of rotation R}
\usage{
  angle(Rs)

  \method{angle}{SO3} (Rs)

  \method{angle}{Q4} (Qs)

  \method{angle}{EA} (eur)
}
\arguments{
  \item{Rs}{rotation matrix, unit quaternion or Euler
  angle}
}
\value{
  angle of rotation

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  Extract angle from rotation.
}
\seealso{
  \code{\link{axis}}
}


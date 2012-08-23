\name{riedist.Q4}
\alias{riedist.Q4}
\title{Riemannian Distance Between Two Random Rotations in quaternion}
\usage{
  riedist.Q4(q, Q = c(1, 0, 0, 0))
}
\arguments{
  \item{q}{The estimate of the central direction}

  \item{Q}{The true central direction}
}
\value{
  the Riemannian distance between the two rotations given
  in quaternion form
}
\description{
  This function will calculate the riemannian distance
  between an estimate of the central direction (in matrix
  or vector form) and the central direction.  By default
  the central direction is taken to be the identity matrix,
  but any matrix in SO(3) will work.  It calls the matrix
  log and matrix exponential functions also given here.
}
\examples{
r<-rvmises(20,0.01)
Qs<-genR(r,space="Q4")
Qp<-mean(Qs)
riedist.Q4(Qp)
}


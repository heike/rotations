\name{genR}
\alias{genR}
\title{Generate Rotations}
\usage{
  genR(r, S = NULL, space = "SO3")
}
\arguments{
  \item{r}{vector of angles}

  \item{S}{The principle direction}

  \item{space}{Indicates the desired representation: matrix
  in SO3, quaternion, or Euler angles}
}
\value{
  a matrix where each row is a sample point in the desired
  space
}
\description{
  Generate rotations according to the Uniform-Axis Random
  Spin methodology
}
\details{
  Given a vector \eqn{u\in\mathbb{R}^2} of length one and
  angle of rotation r, a rotation can be formed using
  Rodrigues formula \deqn{\cos(r)I_{3\times
  3}+\sin(r)\Phi(u)+(1-\cos(r))uu^\top} where \eqn{\Phi(u)}
  is a \eqn{3\times 3} skew-symmetric matirix with upper
  triangular elements \eqn{-u_3}, \eqn{u_2} and \eqn{-u_1}
  in that order.
}
\examples{
r<-rvmises(20,0.01)
genR(r)
}
\references{
  Bingham M, Nordman D and Vardeman S (2009). "Modeling and
  Inference for Measured Crystal Orientations and a
  Tractable Class of Symmetric Distributions for Rotations
  in three Dimensions." _Journal of the American
  Statistical Association_, *104*(488), pp. 1385-1397.
}


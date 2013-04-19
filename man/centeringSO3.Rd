\name{centeringSO3}
\alias{centeringSO3}
\title{Centering function}
\usage{
  centeringSO3(Rs, S)
}
\arguments{
  \item{Rs}{the sample to be centered}

  \item{S}{the rotation to sample around}
}
\value{
  The centered sample
}
\description{
  This function will take the sample Rs and return teh
  sample Rs centered at S, i.e., the returned sample is
  S'Rs, so if S is the true center then the projected mean
  should be id.SO3
}


\name{trim}
\alias{trim}
\title{Function to Trim the Sample}
\usage{
  trim(Rs, alpha)
}
\arguments{
  \item{Rs}{The sample of random rotations}

  \item{alpha}{The percent of observations to be trimeed}
}
\value{
  S3 \code{trim} object; a sample of size n-(n*alpha) of
  random roatations
}
\description{
  This function will take a sample of size n random
  rotations, find the observations beyond the alpha/2
  percentile and delete them from the sample. To determine
  which observations to remove, the average distance of
  each sample point from the estimated central direction is
  used.
}


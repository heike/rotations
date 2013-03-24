\name{pointsXYZ}
\alias{pointsXYZ}
\title{Project rotation data onto sphere}
\usage{
  pointsXYZ(data, center, column = 1)
}
\arguments{
  \item{data}{data frame of rotation matrices in 3 x 3
  matrix representation}

  \item{center}{point about which to center the
  observations}

  \item{column}{integer 1 to 3 indicating which column to
  display}
}
\value{
  data frame with columns X, Y, Z standing for the
  respective coordinates in 3d space
}
\description{
  Projection of rotation matrices onto sphere with given
  center.
}


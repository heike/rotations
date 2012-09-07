\name{eyeBall}
\alias{eyeBall}
\title{A novel approach to visualizing random rotations.}
\usage{
  eyeBall(Rs, center = id.SO3, column = 1,
    show_estimates = FALSE, xlimits = c(-1, 1),
    ylimits = c(-1, 1), ...)
}
\arguments{
  \item{Rs}{the sample of n random rotations}

  \item{center}{point about which to center the
  observations}

  \item{column}{integer 1 to 3 indicating which column to
  display}

  \item{show_estimates}{rather to display the four
  estimates of the principal direction or not}

  \item{xlimits}{limits for the x-axis, appropriate range
  is [-1,1]}

  \item{ylimits}{limits for the y-axis, appropriate range
  is [-1,1]}

  \item{...}{Additional arguments passed to ggplot2}
}
\value{
  a ggplot2 object with the data dispalyed on a blank
  sphere
}
\description{
  This function produces a three-dimensional globe onto
  which the on column of the provided sample is drawn.  The
  data are centered around a provided matrix and the user
  can choose to display this center or not.  Based on
  \code{ggplot2} package by \cite{Wickham (2009)}.
}
\examples{
r<-rvmises(20,1.0)
Rs<-genR(r)
eyeBall(Rs,center=mean(Rs),show_estimates=TRUE,shape=4)
}
\references{
  Wickham H (2009). _ggplot2: elegant graphics for data
  analysis_. Springer New York. ISBN 978-0-387-98140-6,
  <URL: http://had.co.nz/ggplot2/book>.
}


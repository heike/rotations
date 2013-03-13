\name{drilldata}
\alias{drilldata}
\docType{data}
\title{
Data from Rancourt, Rivest and Asselin (2000)
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(drilldata)}
\format{
  A data frame with 720 observations on the following 8 variables.
  \describe{
    \item{\code{Subject}}{Subject number (1-8)}
    \item{\code{Joint}}{Joint number (1=wrist, 2=elbow, 3=shoulder)}
    \item{\code{Position}}{Drilling position (1-6)}
    \item{\code{Replicate}}{Replicate number (1-5)}
    \item{\code{Q1}}{First element of rotation vector}
    \item{\code{Q2}}{Second element of rotation vector}
    \item{\code{Q3}}{Third element of rotation vector}
    \item{\code{Q4}}{Fourth element of rotation vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(drilldata)
## maybe str(drilldata) ; plot(drilldata) ...
}
\keyword{datasets}

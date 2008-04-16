\name{monoxide}
\alias{monoxide}
\docType{data}
\title{ Carbon monoxide levels at an oil refinery in California }
\description{
  The measurements were made close to a California oil refinery in 1990-1993.
}
\usage{data(monoxide)}
\format{
  A data frame with 9 observations on the following 2 variables.
  \describe{
    \item{\code{CO}}{carbon monoxide level (ppm)}
    \item{\code{date}}{date of the measurement}
  }
}
\details{
  Example dataset in Krishnamoorthy and Mathew (2003), p. 116-118. More
  details on the data can be found at DASL.
}
\source{
  Data from the Data and Story Library (DASL)
  \url{http://lib.stat.cmu.edu/DASL/Stories/MeasuringAirPollution.html}
}
\references{
  Krishnamoorthy, K. and Mathew, Th. (2003). Inferences on the means of
  lognormal distributions using generalized \emph{p}-values and generalized
  confidence intervals, Journal of statistical planning and inference,
  115, 103-121.
  
  Olsson, U. (2005). Confidence intervals for the mean of a log-normal
  distribution, Journal of Statistics Education, 13, 1.

  \url{http://amstat.org.publications/jse/v13n1/olsson.html}
}
\examples{
  data(monoxide)
}
\keyword{datasets}

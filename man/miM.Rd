\name{miM}
\alias{miM}
\title{The miniMax criterion baesd on an approximate fill distance measure}
\usage{
miM(D, num = 50)
}
\arguments{
\item{D}{a design matrix, rows are design locations, columns are design factors}

\item{num}{Optional, default is "50". The fineness of the gridded points to divide the design space. Each dimension is evenly divided by num+1 points. Lower this parameter when dimension is high to reduce computing time.}
}
\value{
fill distance with 4 decimals
}
\description{
The miniMax criterion baesd on an approximate fill distance measure
}
\details{
This function calculates the approximate fill distance for the design by using a set of gridded points, the maximum error of the value can be computed.
}
\examples{

#Compute the approximate fill distance of a design D
d=miM(D = cbind(c(0.875,0.375,0.125,0.625),c(0.375,0.125,0.625,0.875)), num = 20)
}

\keyword{ Computer Experiments Criteria}

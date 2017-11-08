\name{md}
\alias{md}
\title{Minimum distance between any two points in the design}
\usage{
md(D)
}
\arguments{
\item{D}{a design matrix, rows are design locations, columns are design factors}
}
\value{
\item{MinimumDistance}{Minimum distance between any two points in the design}

\item{number}{number of pairs in the design achieve the minimum distance}
}
\description{
Minimum distance between any two points in the design
}
\examples{

#compute the minimum distance between any two points in design D
d=md(D = cbind(c(0.875,0.375,0.125,0.625),c(0.375,0.125,0.625,0.875)))
d$MinimumDistance
d$number
}

\keyword{ Computer Experiments Criteria }

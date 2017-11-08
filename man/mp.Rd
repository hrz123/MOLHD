\name{mp}
\alias{mp}
\title{Computer the MaxPro Criterion for a design.}
\usage{
mp(D)
}
\arguments{
\item{D}{a design matrix}
}
\value{
The MaxPro Criterion with 4 decimals
}
\description{
Computer the MaxPro Criterion for a design.
}
\details{
This function is to compute the MaxPro criterion for measuring projection characteristic of a computer experiment.
}
\examples{

#compute the mp criterion of a random LHD(20,2)
D=LHD(20,2)$standDesign
mp(D)
}

\keyword{ Computer Experiments Criteria }

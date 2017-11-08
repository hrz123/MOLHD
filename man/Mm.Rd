\name{Mm}
\alias{Mm}
\title{Computer the approximate Maximin Criterion for a design.}
\usage{
Mm(D, power = 100)
}
\arguments{
\item{D}{a design matrix}

\item{power}{Optional, default is "100". The power parameter r in the average reciprocal inter-point distance measure. When r is approaching infinity, minimizing the average reciprocal inter-point distance measure is equivalent to maximizing the minimum distance among the design points.}
}
\value{
The approximate Maximin criterion with 4 decimals
}
\description{
Computer the approximate Maximin Criterion for a design.
}
\examples{

#Compute the maximin criterion of a random LHD(20,2)
des=LHD(n = 20,p = 2)$standDesign
Mm(D=des, power=150)
}

\keyword{ Computer Experiments Criteria }

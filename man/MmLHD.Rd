\name{MmLHD}
\alias{MmLHD}
\title{Generate the optimal Maximin Latin Hypercube Design.}
\usage{
MmLHD(n, p, power = 100, temp0 = 0, nstarts = 1, times = 300,
  maxiter = 1e+06)
}
\arguments{
\item{n}{number of runs desired}

\item{p}{number of variables desired}

\item{power}{Optional, default is "100". The power parameter r in the average reciprocal inter-point distance measure. When r turns to infinity, minimizing the average reciprocal inter-point distance measure is equivalent to maximizing the minimum distance among the design points.}

\item{temp0}{Initial temperature}

\item{nstarts}{Optional, default is "1". The number of random starts}

\item{times}{Optional, default is "300". The maximum number of non-improving searches allowed. Lower this parameter if you expect the search to converge faster.}

\item{maxiter}{Optional, default is "1e+06".The maximum total number of iterations for each random start. Lower this number if the design is prohibitively large and you want to terminate the algorithm prematurely to report the best design found}
}
\value{
\item{design}{The optimal Maximin design matrix}

\item{criterion}{The opproximate Maximin criterion of the design under chosen "power" parameter}

\item{iterations}{The total iterations}

\item{time_rec}{Time to complete the search}
}
\description{
Generate the optimal Maximin Latin Hypercube Design.
}
\details{
This function is to search the optimal Maximin design using columnwise exchange algorithm coupled with the simulated annealing algorithm and several computational shortcuts to improve efficiency.
}
\examples{

#Generate the optimal maximin distance LHD(20,2)
D=MmLHD(n=20,p=2)
D$design
D$criterion
}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }

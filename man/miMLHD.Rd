\name{miMLHD}
\alias{miMLHD}
\title{Generate the optimal Latin Hypercube Design based on the miniMax criterion.}
\usage{
miMLHD(n, p, num = 50, temp0 = 0, nstarts = 1, times = 300,
  maxiter = 1e+06)
}
\arguments{
\item{n}{number of runs desired}

\item{p}{number of variables desired}

\item{num}{Optional, default is "50". The fineness of the gridded points to divide the design space. Each dimension is evenly divided by num+1 points. Lower this parameter when dimension is high to reduce computing time.}

\item{temp0}{Initial temperature for simulated annealing}

\item{nstarts}{Optional, default is "1". The number of random starts}

\item{times}{Optional, default is "300". The maximum number of non-improving searches allowed before terminating the search.}

\item{maxiter}{Optional, default is "1e+06".The maximum total number of iterations for each random start. Lower this number if the design is prohibitively large and you want to terminate the algorithm prematurely to report the best design found}
}
\value{
\item{design}{The optimal miniMax design matrix}

\item{criterion}{The opproximate miniMax criterion for the chosen fineness of the grids}

\item{iterations}{The total iterations}

\item{time_rec}{Time to complete the search}
}
\description{
Generate the optimal Latin Hypercube Design based on the miniMax criterion.
}
\details{
This function is to search the optimal Latin Hypercube design based on the miniMax criterion using the columnwise exchange algorithm coupled with the simulated annealing algorithm, and several computational shortcuts to improve efficiency. The approximate miniMax criterion is computed by using a set of gridded points to approximate the continuous design space, the maximum error of the value can be computed.(Can only work in relatively low dimensions)
}
\examples{

#Generate the optimal minimax distance LHD(20,2)
\dontrun{
D=miMLHD(n=20,p=2)
D$design
D$criterion
}
}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }

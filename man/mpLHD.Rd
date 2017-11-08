\name{mpLHD}
\alias{mpLHD}
\title{Generate the optimal MaxPro Latin Hypercube Design.}
\usage{
mpLHD(n, p, temp0 = 0, nstarts = 1, times = 300, maxiter = 1e+06)
}
\arguments{
\item{n}{number of runs desired}

\item{p}{number of design factors desired}

\item{temp0}{Initial temperature for simulated annealing}

\item{nstarts}{Optional, default is "1". The number of random starts}

\item{times}{Optional, default is "300". The maximum number of non-improving searches allowed. Lower this parameter if you expect the search to converge faster.}

\item{maxiter}{Optional, default is "1e+06".The maximum total number of iterations. Lower this number if the design is prohibitively large and you want to terminate the search prematurely to report the best design found}
}
\value{
\item{design}{The optimal LHD design matrix based on the MaxPro criterion}

\item{criterion}{The MaxPro criterion of the selected optimal LHD design}

\item{iterations}{The total iterations}

\item{time_rec}{Time to complete the search}
}
\description{
Generate the optimal MaxPro Latin Hypercube Design.
}
\details{
This function is to search the optimal Latin Hypercube Design based on the MaxPro criterion using the columnwise exchange algorithm coupled with the simulated annealing algorithm, and several computational shortcuts to improve efficiency.
}
\examples{

#Generate a optimal maximum projection LHD(20,2) design
D=mpLHD(n=20,p=2)
D$design
D$criterion
}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }

\name{pfMm}
\alias{pfMm}
\title{Generate the Pareto front for the optimal Latin Hypercube Designs based on both the Maximin and miniMax criteria.}
\usage{
pfMm(n, p, crlim, num, nstarts = 1, times = 300, maxiter = 1e+06,
  temp0 = 0, wtset = cbind(c(1, 0), c(0.8, 0.2), c(0.6, 0.4), c(0.4, 0.6),
  c(0.2, 0.8), c(0, 1)))
}
\arguments{
\item{n}{number of runs desired}

\item{p}{number of design factors desired}

\item{crlim}{a matrix saving the best and worst values for each criterion to be used for defining the scaling choices for converting the natural criteria values onto a desirability scale between 0 and 1. Each column corresponds to one criterion. The best and worst values are recommended based on the values from each single criterion search. It is recommended that slightly wider range is used for defining the scaling choice for the Pareto front search.}

\item{num}{The fineness of the grids to approximiate the approximate the continuous design space. Lower this parameter when dimension is high to reduce computing time.}

\item{nstarts}{Optional, default is "1". The number of random starts}

\item{times}{Optional, default is "300". The maximum number of non-improving searches allowed before terminating the search. Lower this parameter if you expect the search to converge faster.}

\item{maxiter}{Optional, default is "1e+06".The maximum total number of iterations. Lower this number if the design is prohibitively large.}

\item{temp0}{Initial temperature for simualted annealing}

\item{wtset}{Optional, default is "cbind(c(1,0),c(0.8,0.2),c(0.6,0.4),c(0.4,0.6),c(0.2,0.8),c(0,1))". The set of weight combinations to guide the search in varied directions. Each column is a weight vector that guides the search in a certain direction.}
}
\value{
\item{pfdes}{The column bind of Pareto front designs whose criteria values are on the Pareto front.}

\item{pfvals}{The Pareto front of criteria values based on the Maximin and miniMax criteria. Columns are the optimization criteria.}

\item{time_rec}{Time to complete the search}
}
\description{
Generate the Pareto front for the optimal Latin Hypercube Designs based on both the Maximin and miniMax criteria.
}
\details{
This function is to search for the Pareto front and the Pareto set of LHDs based on the Maximin and miniMax criteria. Each design on Pareto front is not dominated by any other design.This function utilizes a version of simulated annealing algorithm and several computational shortcuts to efficiently generate the optimal Latin hypercube designs. Choose a lower maximum limit of the criteria but high enough for Pareto front designs will save the computing time.
}
\examples{

#Generate the Pareto front designs of maximin and minimax distance criterion for LHD(10,2)
\dontrun{
D1=MmLHD(n=10,p=2,nstarts=30)
D2=miMLHD(n=10,p=2,num=15,nstarts=30)
Mmlim=c(D1$criterion-0.2,D1$criterion-0.2+2)
mMlim=c(D2$criterion-0.05,D2$criterion-0.05+0.5)
crlim=cbind(Mmlim,mMlim)
pf=pfMm(10,2,crlim,num = 15,nstarts = 30)
pf$pfvals
pf$pfdes
pf$time_rec
}
}
\keyword{ Design of Experiments }
\keyword{ Computer Experiments }

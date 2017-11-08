\name{cpf2}
\alias{cpf2}
\title{Combine Pareto front designs of 2 criteria}
\usage{
cpf2(newdes, newpfval, curdes, curpfval)
}
\arguments{
\item{newdes}{a matrix which is a column bind of new designs}

\item{newpfval}{a matrix each row is 2 criteria correponding to each design}

\item{curdes}{a matrix which is a column bind of current designs on Pareto front}

\item{curpfval}{a matrix each row is 2 criteria correponding to each Pareto front design}
}
\value{
\item{pfdes}{The column bind of Pareto front designs}

\item{pfvals}{The Pareto front values corresponding to the Pareto front designs}
}
\description{
Combine Pareto front designs of 2 criteria
}
\details{
This function is used to combine 2 criteria Pareto front designs
}
\examples{

#Combine Pareto fronts each with 5 random starts for Mm and mp criteria
\dontrun{
pf1=pfMp(20,2,crlim = cbind(c(4.5,6.5),c(26,36)),nstarts = 5)
pf2=pfMp(20,2,crlim = cbind(c(4.5,6.5),c(26,36)),nstarts = 5)
pfnew=cpf2(pf1$pfdes,pf1$pfvals,pf2$pfdes,pf2$pfvals)
pfnew$pfdes
pfnew$pfvals
}
}
\keyword{ Pareto front }

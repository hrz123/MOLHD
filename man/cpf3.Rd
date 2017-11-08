\name{cpf3}
\alias{cpf3}
\title{Combine Pareto front designs of 3 criteria}
\usage{
cpf3(newdes, newpfval, curdes, curpfval)
}
\arguments{
\item{newdes}{a matrix which is a column bind of new designs}

\item{newpfval}{a matrix each row is 3 criteria correponding to each design}

\item{curdes}{a matrix which is a column bind of current designs on Pareto front}

\item{curpfval}{a matrix each row is 3 criteria correponding to each Pareto front design}
}
\value{
\item{pfdes}{The column bind of Pareto front designs}

\item{pfvals}{The Pareto front values corresponding to the Pareto front designs}
}
\description{
Combine Pareto front designs of 3 criteria
}
\details{
This function is used to combine 3 criteria Pareto front designs
}
\examples{

#Combine Pareto fronts each with 1 random start for Mm, mp and mM criteria
\dontrun{
pf1=pfMpm(20,2,crlim = cbind(c(4.5,6.5),c(26,36),c(0.12,0.62)),num = 15,nstarts = 1)
pf2=pfMpm(20,2,crlim = cbind(c(4.5,6.5),c(26,36),c(0.12,0.62)),num = 15,nstarts = 1)
pfnew=cpf3(pf1$pfdes,pf1$pfvals,pf2$pfdes,pf2$pfvals)
pfnew$pfdes
pfnew$pfvals
}
}

\keyword{ Pareto front }

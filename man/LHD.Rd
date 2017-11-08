\name{LHD}
\alias{LHD}

\title{Generate a random Latin Hypercube design}
\usage{
LHD(n, p)
}
\arguments{
\item{n}{number of runs desired}

\item{p}{number of design factors}
}
\value{
\item{design}{a Latin Hypercube Design that is not scaled (i.e. the grid point locations are integers)}

\item{standDesign}{a standard Latin Hypercube Design that is scaled to (0,1); design locaitons are placed at the centers of selected grids.}
}
\description{
Generate a random Latin Hypercube design
}
\examples{

#Generate a random Latin hypercube design with 20 runs and 2 variables
D<-LHD(n = 20,p = 2)
D$design
D$standDesign
}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }

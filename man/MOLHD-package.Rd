\name{MOLHD-package}
\alias{MOLHD-package}
\alias{MOLHD}
\docType{package}

\title{Multiple objective Latin hypercube design}
\description{
The MOLHD package provides useful and efficient functions for generating the optimal Maximin distance, Maximum Projection and miniMax distance (only for low dimensions) designs within the class of Latin hypercube designs for computer experiments. Ant it provides functions generating Pareto front optimal designs for each two of the three criteria and all the three criteria within the class of Latin hypercube designs. It also provides functions to compute the criteria for a given design.
}
\details{
\tabular{ll}{
Package: \tab MOLHD\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2017-11-06\cr
Lisense: \tab LGPL\cr
}
This package contains functions for generating the optimal maximin distance designs, maximum projection designs and minimax distance designs for low dimensions within the class of Latin hypercube designs (LHDs). This packages also contains functions for generating designs on the Pareto front of each two of the three criteria as maximin distance criterion, minimax distance criterion, and maximum projection criterion. This package also contains functions to compute each criterion for a random Latin hypercube design.

Since minimax distance design is computational expensive, it is only approximately extimated when the design is at low dimension.
}

\author{
Ruizhe Hou, Lu Lu

Maintainer: Ruizhe Hou<houruizhe1210@gmail.com>
}
\references{
Morris, M. D. and Mitchell, T. J. (1995), "Exploratory Designs for Computation Experiments," \emph{Journal of Statistical Planning and Inference}. <doi:10.1016/0378-3758(94)00035-T>

Lu Lu and Christine M. Anderson-CookTimothy J. Robinson (2011), "Optimization of Designed Experiments Based on Multiple Criteria Utilizing a Pareto Frontier," \emph{Technometrics}.
<doi:10.1198/Tech.2011.10087>

Joseph, V. R., Gul, E., and Ba, S. (2015), "Maximum Projection Designs for Computer experiments," \emph{Biometrika}. <doi:10.1093/biomet/asv002>

}
\keyword{ package }

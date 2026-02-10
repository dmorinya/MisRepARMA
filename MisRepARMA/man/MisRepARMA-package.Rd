\encoding{utf8} 
\name{MisRepARMA-package}
\alias{MisRepARMA-package}
\alias{MisRepARMA}
\docType{package}
\title{
Misreported time series analysis
}
\description{
Provides a simple and trustworthy methodology for the analysis of misreported continuous time series. See Moriña, D, Fernández-Fontelo, A, Cabaña, A, Puig P. (2021) <https://arxiv.org/abs/2003.09202v2>.
}
\details{
\tabular{ll}{
Package: \tab MisRepARMA\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2021-07-14\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
The package implements function \code{\link{fitMisRepARMA}}, 
which is able to fit an ARMA time series model to misreported data, and the function 
\code{\link{reconstruct}} which is able to reconstruct the most likely real series.
}

\author{
David Moriña, Amanda Fernández-Fontelo, Alejandra Cabaña, Pedro Puig

Mantainer: David Moriña Soler <dmorina@ub.edu>
}

\seealso{
\code{\link{MisRepARMA-package}}, \code{\link{fitMisRepARMA}}, \code{\link{reconstruct}}
}

\references{
Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Kunsch, H.R. (1989) The jackknife and the bootstrap for general stationary observations. Annals of Statistics, \bold{17}, 1217–1241.

Moriña, D., Fernández-Fontelo, A., Cabaña, A., Puig, P. (2021): New statistical model for misreported data with application to current public health challenges. arXiv preprint (https://arxiv.org/pdf/2003.09202.pdf)

Politis, D.N. and Romano, J.P. (1994) The stationary bootstrap. Journal of the American Statistical Association, \bold{89}, 1303–1313.
}
\keyword{ MisRepARMA }


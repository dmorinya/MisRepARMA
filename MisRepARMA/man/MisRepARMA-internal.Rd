\encoding{utf8} 
\name{MisRepARMA-internal}
\alias{MisRepARMA-internal}
\alias{estimate}
\alias{ran.genf}
\alias{summary.fitMisRepARMA}
\alias{print.summary.fitMisRepARMA}

\title{Internal MisRepARMA functions}
\description{Internal MisRepARMA functions}
\usage{
estimate(data, tol, p_AR, q_MA, covars=NULL, misReport="U")
ran.genf(data, n, ran.args)
\method{summary}{fitMisRepARMA}(object, ...)
\method{print.summary}{fitMisRepARMA}(x, ...)
}

\details{These functions are not to be called by the user}
\seealso{
\code{\link{MisRepARMA-package}}, \code{\link{fitMisRepARMA}}, \code{\link{reconstruct}}
}

\keyword{ internal }

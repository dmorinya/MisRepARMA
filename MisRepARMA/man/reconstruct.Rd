\encoding{utf8}
\name{reconstruct}
\alias{reconstruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reconstruct the most likely series
}
\description{
  Reconstructs the most likely series.
}
\usage{
  reconstruct(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  object of class \code{fitMisRepARMA}.
}
}

\value{
  the function returns a vector of the same length of \code{data} containing the reconstruction of the most likely series.
}

\references{
D. Moriña, A. Fernández-Fontelo, A. Cabaña, P. Puig (2021): New statistical model for misreported data with application to current public health challenges. arXiv preprint (https://arxiv.org/pdf/2003.09202.pdf)

Davison, A. C. and Hinkley, D. V. (1997) Bootstrap Methods and Their Applications. Cambridge University Press, Cambridge. ISBN
  0-521-57391-2
}

\author{
David Moriña, Amanda Fernández-Fontelo, Alejandra Cabaña, Pedro Puig
}

\seealso{
\code{\link{MisRepARMA-package}}, \code{\link{fitMisRepARMA}}
}

\examples{
### Simulate underreported time series data
x <- arima.sim(model=list(ar=0.4), n=50)
ind <- rbinom(50, 1, 0.6)
y <- ifelse(ind==0, x, x*0.3)
pr <- fitMisRepARMA(y, 1e-8, 5, 0.05, 1, 0, covars=NULL, misReport="U", method="freq")
x <- reconstruct(pr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MisRepARMA }

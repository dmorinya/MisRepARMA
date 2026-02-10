\encoding{utf8}
\name{fitMisRepARMA}
\alias{fitMisRepARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit ARMA model to misreported time series data
}
\description{
  Fits an ARMA model to misreported time series data.
}
\usage{
  fitMisRepARMA(y, tol, B, p_AR, q_MA, covars=NULL, misReport="U", method="freq", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  a numeric vector or time series giving the original data.
}
  \item{tol}{
  tolerance limit to stop the iterative algorithm.
}
  \item{B}{
  the number of bootstrap series to compute.
}
\item{p_AR}{
  order of the AR part.
}
\item{q_MA}{
  order of the MA part.
}
\item{covars}{
  matrix of explanatory variables. Its default value is \code{NULL}.
}
\item{misReport}{
  direction of misreporting issue. Its default value is \code{U} for underreported data, can also take the value \code{O} for overreported data.
}
\item{method}{
  method to estimate the parameters. Its default value is \code{freq} for frequentist approach, can also take the value \code{bayes} for Bayesian approach.
}
\item{...}{
  additional arguments to pass to \code{tsboot}, for instance those regarding parallelization.
}
}

\details{
The model based resampling scheme with \code{B} bootstrap resamples is computed. This
}

\value{
An object of class \code{fitMisRepARMA} with the following elements is returned:

\itemize{
  \item \code{data}: The original time series.
  \item \code{t0}: The results of applying statistic to the original series.
  \item \code{t}: Estimates on each replicated time series.
  \item \code{call}: The original call to tsboot.
}
}

\author{
David Moriña, Amanda Fernández-Fontelo, Alejandra Cabaña, Pedro Puig
}

\references{
Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Kunsch, H.R. (1989) The jackknife and the bootstrap for general stationary observations. Annals of Statistics, \bold{17}, 1217–1241.

Moriña, D., Fernández-Fontelo, A., Cabaña, A., Puig, P. (2021): New statistical model for misreported data with application to current public health challenges. arXiv preprint (https://arxiv.org/pdf/2003.09202.pdf)

Politis, D.N. and Romano, J.P. (1994) The stationary bootstrap. Journal of the American Statistical Association, \bold{89}, 1303–1313.
}

\seealso{
\code{\link{MisRepARMA-package}}, \code{\link{reconstruct}}
}

\examples{
### Simulate underreported time series data
set.seed(12345)
x <- arima.sim(model=list(ar=0.4), n=50)
ind <- rbinom(50, 1, 0.6)
y <- ifelse(ind==0, x, x*0.3)
mod <- fitMisRepARMA(y, 1e-6, 3, 0.05, 1, 0, covars=NULL, misReport="U", method="freq")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MisRepARMA }

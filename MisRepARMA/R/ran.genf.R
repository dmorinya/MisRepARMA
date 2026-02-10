ran.genf <- function(data, n, ran.args)
{
  if (length(ran.args$ar)==0)
  {
    x <- as.numeric(arima.sim(model=list(ma=ran.args$ma), n=n, 
                   innov=rnorm(n, ran.args$intercept, sqrt(ran.args$var))))
  }
  if (length(ran.args$ma)==0)
  {
    x <- as.numeric(arima.sim(model=list(ar=ran.args$ar), n=n, 
                   innov=rnorm(n, ran.args$intercept, sqrt(ran.args$var))))
  }
  if (length(ran.args$ar)!=0 & length(ran.args$ma)!=0)
  {
    x <- as.numeric(arima.sim(model=list(ar=ran.args$ar, ma=ran.args$ma), n=n, 
                              innov=rnorm(n, ran.args$intercept, sqrt(ran.args$var))))
  }
  y <- ifelse(ran.args$z==0, x, x*ran.args$q)
  return(y)
}

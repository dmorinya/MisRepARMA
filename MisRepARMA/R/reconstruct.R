reconstruct <- function(object)
{
  if (class(object)!="fitMisRepARMA") stop("Reconstruction should be done on a fitMisRepARMA object.")
  aux   <- attr(object, "z")
  q <- attr(object, "q")
  if (is.null(attr(object, "covars"))) 
  {
    x_rec <- ifelse(aux==0, object$data, object$data/q)
  }else{
    data2 <- object$data - attr(object, "covars")
    x_rec <- ifelse(aux==0, data2, data2/q) + attr(object, "covars")
  }
  return(x_rec)
}

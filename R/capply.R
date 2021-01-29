##' Apply function over tables
##' 
##' Apply a function to each contingency table in a \code{tables} object.
##' 
##' @param x object of class \code{tables}
##' @param f function to apply to each table
##' @param ... additional arguments to \code{f}
##' 
##' @return a vector, matrix or list of outputs from the function \code{f}.
##' 
##' @export capply
capply <- function(x, f, ...) {
  tmp <- apply(as.array(x), 1, f, ...)
  if (is.matrix(tmp)) return(t(tmp))
  else tmp
}

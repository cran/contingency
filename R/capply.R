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
##' @export
capply <- function(x, f, ...) {
  tmp <- apply(as.array(x), 1, f, ...)
  if (is.matrix(tmp)) return(t(tmp))
  else tmp
}

##' Turn distributions into tables
##' 
##' @param n number of distributions to generate
##' @param f function that generates a probability distribution
##' @param ... arguments to \code{f}
##' 
##' @return a tables object containing the outputs of \code{f}
##' 
##' @export
repTables <- function (n, f, ...) {
  args <- list(...)
  out <- replicate(n, do.call(f, args))
  if (is.null(dim(out))) stop("f should give distributions")
  
  out <- aperm(out, c(length(dim(out)), seq_along(dim(out)[-1])))
  out <- as_tables(out, tdim=dim(out)[-1])
  
  out
}

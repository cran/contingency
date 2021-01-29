##' Check conditional independence
##' 
##' Gives a numerical check that a (conditional) independence holds in a probability distribution.
##' 
##' @details just tests to an appropriate numerical precision that a conditional independence
##' holds: this is *not* a statistical test for conditional independence.
##' If \code{A} and \code{B} overlap with \code{C} then these vertices are ignored.  If \code{A}
##' and \code{B} intersect with one another (but not \code{C}) then the solution is always
##' false.
##' 
##' @param x an array or object of class \code{tables}
##' @param A,B the sets of variables whose independence is to be tested
##' @param C conditioning set (possibly empty)
##' @param eps tolerance parameter
##' @param ... other arguments to methods
##' 
##' @return A logical, or a vector of logicals of the same length as the number
##' of tables provided, indicating whether the conditional independence seems 
##' to hold numerically.
##' 
##' @export checkCI
checkCI <- function(x, A, B, C=integer(0), eps=.Machine$double.eps, ...) {
  UseMethod("checkCI")
}

#' @describeIn checkCI method for \code{array} object
#' @method checkCI array
#' @export 
checkCI.array <- function(x, A, B, C=integer(0), eps=.Machine$double.eps, ...) {
  A <- setdiff(A,C)
  B <- setdiff(B,C)
  if (length(A) == 0 || length(B) == 0) return(TRUE)

  dA <- prod(dim(x)[A])
  dB <- prod(dim(x)[B])
  dC <- prod(dim(x)[C])
  
  x <- margin(x, c(C,A,B))
  dim(x) <- c(dC, dA, dB)
  
  tmp <- apply(x, 1, function(x) svd(x, nv=0, nu=0)$d[2])
  if (any(tmp > eps)) return(FALSE)
  
  return(TRUE)
}

#' @describeIn checkCI method for \code{tables} object
#' @method checkCI tables
#' @export
checkCI.tables <- function(x, A, B, C=integer(0), eps=.Machine$double.eps, ...) {
  n <- ntables(x)
  A <- setdiff(A,C)
  B <- setdiff(B,C)
  if (length(A) == 0 || length(B) == 0) return(rep(TRUE, n))
  if (length(intersect(A,B)) > 0) return(rep(FALSE, n))
  
  dA <- prod(tdim(x)[A])
  dB <- prod(tdim(x)[B])
  dC <- prod(tdim(x)[C])
  
  x <- margin(x, c(C,A,B))
  tdim(x) <- c(dC, dA, dB)
  
  out <- logical(n)
  
  for (i in seq_len(n)) {
    tmp <- apply(x[i,], 1, function(x) svd(x, nv=0, nu=0)$d[2])
    out[i] <- all(tmp <= eps)
  }
  
  return(out)
}

# ##' Enforce Conditional Independence
# ##' 
# ##' Returns the MLE (i.e. KL minimizer) for a model defined by a 
# ##' single conditional independence.
# ##' 
# ##' @param x
# ##' @param A,B,C margins of \code{x}
# ##' 
# ##' @details Given p(A, B, C, D) this returns p(A, C)*p(B | C)*p(D | A, B, C).
# ##' 
# fitCI <- function(x, A, B, C=integer(0), ...) {
#   UseMethod("fitCI")
# }
# 
# fitCI.default <- function(x, A, B, C=integer(0), ...) {
#   A <- setdiff(A,C)
#   B <- setdiff(B,C)
#   if (length(A) == 0 || length(B) == 0) return(x)
#   if (length(intersect(A, B)) > 0) stop("Separator sets A and B overlap")
#   
#   D = setdiff(seq_along(dim(x)), c(A,B,C))
#   margin2(x, c(A,C))*conditional2(x, B, C)*conditional2(x, D, c(A,B,C))
# }
# 
# fitCI.tables <- function(x, A, B, C=integer(0), ...) {
#   A <- setdiff(A,C)
#   B <- setdiff(B,C)
#   if (length(A) == 0 || length(B) == 0) return(x)
#   if (length(intersect(A, B)) > 0) stop("Separator sets A and B overlap")
#   
#   D = setdiff(seq_along(tdim(x)), c(A,B,C))
#   margin2(x, c(A,C))*conditional2(x, B, C)*conditional2(x, D, c(A,B,C))
# }

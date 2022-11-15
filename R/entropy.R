##' Calculate entropy of discrete distribution
##' 
##' @param p non-negative numeric vector
##' @param ... other arguments to methods
##' 
##' @return A numeric value of the entopy, or vector of entropies.
##' 
##' @export entropy
entropy <- function(p, ...) {
  UseMethod("entropy")
}

##' @describeIn entropy Default method for vectors
##' @method entropy default
##' @export
entropy.default <- function(p, ...) {
  q <- p[p>0]
  -sum(q*log(q))/sum(q)
}

##' @describeIn entropy Method for arrays
##' @method entropy array
##' @export
entropy.array <- function(p, margin, ...) {
  if (!missing(margin)) p <- marginTable(p, margin)
  q <- p[p > 0]
  -sum(q*log(q))/sum(q)
}

##' @describeIn entropy Method for \code{tables} object
##' @param margin margin to consider
##' @method entropy tables
##' @export
entropy.tables <- function(p, margin, ...) {
  if (!missing(margin)) p <- margin.tables(p, margin)
  tmp <- p*log(p)
  tmp[is.nan(tmp)] <- 0
  -rowSums(tmp)/rowSums(p)
}

##' (Conditional) mutual information
##' 
##' @param p numeric array or \code{tables} class
##' @param m1,m2 margins for mutual information
##' @param condition conditional margin
##' @param ... other arguments to methods
##' 
##' @return Numeric value for mutual information, or a vector of mutual 
##' information values.
##' 
##' @export mutualInf
mutualInf <- function(p, m1, m2, condition, ...) {
  UseMethod("mutualInf")
}

##' @describeIn mutualInf Default method for vectors
##' @method mutualInf default
##' @export
mutualInf.default <- function(p, m1, m2, condition, ...) {
  if (missing(condition)) condition <- integer(0)
  if (length(intersect(m1,m2)) > 0 || length(intersect(c(m1,m2),condition))) stop("Variable sets must be disjoint")
  else tmp <- margin(p, c(m1, m2, condition))/sum(p)
  
  d1 <- prod(dim(tmp)[seq_along(m1)])
  d2 <- prod(dim(tmp)[seq_along(m2)+length(m1)])
  d3 <- prod(dim(tmp))/(d1*d2)
  dim(tmp) <- c(d1,d2,d3)

  p13 <- conditionTable2(tmp, c(1,3), c())
  p23 <- conditionTable2(tmp, c(2,3), c())
  p3 <- conditionTable2(tmp, c(3), c())

  out <- tmp*(log(tmp) - log(p13) - log(p23) + log(p3))
  return(sum(out[tmp > 0]))
}

##' @describeIn mutualInf Method for \code{tables} object
##' @method mutualInf tables
##' @export
mutualInf.tables <- function(p, m1, m2, condition, ...) {
  if (missing(condition)) condition <- integer(0)
  if (length(intersect(m1,m2)) > 0 || length(intersect(c(m1,m2),condition))) stop("Variable sets must be disjoint")
  else tmp <- margin(p, c(m1, m2, condition))
  tmp <- tmp/rowSums(tmp)
  
  ds <- tdim(tmp)
  d1 <- prod(ds[seq_along(m1)])
  d2 <- prod(ds[seq_along(m2)+length(m1)])
  d3 <- prod(ds)/(d1*d2)
  tdim(tmp) <- c(d1,d2,d3)
  
  p13 <- conditional2(tmp, c(1,3), c())
  p23 <- conditional2(tmp, c(2,3), c())
  p3 <- conditional2(tmp, c(3), c())
  
  out <- tmp*(log(tmp) - log(p13) - log(p23) + log(p3))
  out[is.nan(out)] = 0
  return(rowSums(out))
}

##' Interaction information
##' 
##' @param p object to find interaction information for
##' @param ... other arguments to methods
##' 
##' @return Numeric value for interaction information, or a vector of 
##' interaction information values.
##' 
##' @export interactionInf
interactionInf <- function(p, ...) UseMethod("interactionInf")

##' @describeIn interactionInf Default method for vectors
##' @param condition variables on which to condition
##' @method interactionInf default
##' @export
interactionInf.default <- function(p, ..., condition) {
  dots <- list(...)
  if (missing(condition)) condition=integer(0)
  if (length(dots) == 1) {
    return(-entropy(p, c(dots[[1]], condition)) + entropy(p, condition))
  }
  else {
    args1 <- c(list(p=p), dots[-1], list(condition=condition))
    args2 <- c(list(p=p), dots[-1], list(condition=c(dots[[1]], condition)))
    return(do.call(Recall, args2) - do.call(Recall, args1))
  }
}

##' Kullback-Leibler Divergence
##' 
##' Get the KL Divergence between two discrete distributions
##' 
##' @param x,y vectors (of probabilities)
##' @param ... other arguments to methods
##' 
##' @return a numberic value, vector or matrix of KL-divergences.
##' 
##' @export kl
kl <- function(x, y, ...) {
  UseMethod("kl")
}

##' @describeIn kl Default method for vectors
##' @method kl default
##' @export
kl.default <- function(x, y, ...) {
  sum(x[y > 0]*log(x[y > 0]/y[y > 0]))
}

##' @describeIn kl Method for \code{tables} object
##' @method kl tables
##' @export 
kl.tables <- function(x, y, ...) {
  
  if (!("tables" %in% class(y))) y <- as_tables(y, tdim=tdim(x)) 
  
  if ((length(tdim(x)) != length(tdim(y))) || 
      (!all(tdim(x) == tdim(y)))) stop("table dimensions not equal")

  nx = ntables(x)
  ny = ntables(y)
  if (nx < ny) {
    if (ny %% nx != 0) warning("longer object length is not a multiple of shorter object length")
    x = x[rep.int(seq_len(nx), ny),,drop=FALSE]
  }
  else if (nx > ny) {
    if (nx %% ny != 0) warning("longer object length is not a multiple of shorter object length")
    y = y[rep.int(seq_len(ny), nx),,drop=FALSE]
  }
  
  tmp <- x*log(x/y)
  tmp[is.nan(tmp)] = 0
  rowSums(tmp)
}


##' Multiinformation
##' 
##' Get the multiinformation for a discrete distribution
##' 
##' @param x vectors (of probabilities)
##' @param ... other arguments to methods
##' 
##' @return a numberic value, vector or matrix of required multiinformation.
##' 
##' @export
multiInf <- function(x, ...) {
  UseMethod("multiInf")
}

##' @describeIn multiInf Default method for vectors and arrays
##' @param margin margin to find multiinformation for
##' @method multiInf default
##' @export
multiInf.default <- function(x, margin=NULL, ...) {
  
  if (!is.null(margin)) {
    x <- margin(x, margin)
  }
  
  if (is.null(dim(x))) return(0)
  
  out <- x
  out[] <- 1
  
  for (i in seq_along(dim(x))) {
    out <- out*margin2(x, i)
  }
  
  entropy(x/out)
}

##' @describeIn multiInf Method for \code{tables} object
##' @method multiInf tables
##' @export 
multiInf.tables <- function(x, margin=NULL, ...) {
  
  if (!is.null(margin)) x <- margin(x, margin)
  
  out <- x
  out[] <- 1
  
  for (i in seq_along(tdim(x))) {
    out <- out*margin2(x, i)
  }
  
  entropy(x/out)
}


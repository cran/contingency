##' Get margin of a table or tables
##' 
##' 
##' @param x a contingency table or \code{tables} object
##' @param ... a contingency table or \code{tables} object
##' 
##' @details \code{margin2} keeps all dimensions, and 
##' hence results will sum to the number of cells summed over.
##' 
##' @return an object of the same class as \code{x}.  The resulting
##' array, or collection of tables, will contain a marginal, conditional
##' or interventional distribution.
##' 
##' @export 
margin <- function(x, ...) UseMethod("margin")
##' @describeIn margin keep all dimensions
##' @export 
margin2 <- function(x, ...) UseMethod("margin2")
##' @describeIn margin conditional distributions
##' @export 
conditional <- function(x, ...) UseMethod("conditional")
##' @describeIn margin  conditional distributions with all dimensions kept
##' @export 
conditional2 <- function(x, ...) UseMethod("conditional2")
##' @describeIn margin interventional distributions
##' @export 
intervention <- function(x, ...) UseMethod("intervention")

##' @method margin default
##' @export
margin.default <- function(x, margin=NULL, order=TRUE, ...) {
  marginTable(x, margin, order)
}

##' @method margin2 default
##' @export
margin2.default <- function(x, margin=NULL, ...) {
  out <- marginTable(x, margin)
  out <- out[patternRepeat0(margin, dim(x), keep.order=TRUE)]
  dim(out) <- dim(x)
  out
}

##' @method margin2 tables
##' @export
margin2.tables <- function(x, margin=NULL, ...) {
  out <- margin(x, margin)
  out <- out[patternRepeat0(c(1,margin+1), c(ntables(x), tdim(x)), keep.order=TRUE)]
  dim(out) <- dim(x)
  tdim(out) <- tdim(x)
  class(out) <- "tables"
  out  
}

##' @method conditional default
##' @export
conditional.default <- function(x, variables, condition = NULL, condition.value = NULL, undef=NaN, ...) {
  out = conditionTable(x, variables, condition, condition.value, undef=undef)
#  if (!is.nan(undef)) out[is.nan(out)] = undef
  out
}

##' @method conditional2 default
##' @export
conditional2.default <- function(x, variables, condition = NULL, undef=NaN, ...) {
  out = conditionTable2(x, variables, condition, undef=undef)
  if (!is.nan(undef)) out[is.nan(out)] = undef
  out
  
}

##' Get the marginal distributions
##' 
##' @param x an object of class \code{tables}
##' @param margin integer vector giving margin to be calculated (1 for rows, etc.)
##' @param order logical indicating whether resulting indices
##' should be in the same order as stated in \code{margin} 
##' @param ... other arguments to function
##' 
##' @details Calculates marginal distributions for each entry in a \code{probMat}.
##' 
##' @return An object of class \code{tables} consisting of the required marginal 
##' distribution.
##' 
##' @method margin tables
##' @export
margin.tables <- function(x, margin=NULL, order=TRUE, ...) {
  if (!order) margin <- sort.int(margin)

  out <- marginTable(as.array(x), c(1,margin+1), order=order)
  dim(out) <- c(dim(out)[1], length(out)/dim(out)[1])
  attr(out, "tdim") <- tdim(x)[margin]
  class(out) <- "tables"
  
  out  
}

##' @method conditional tables
##' @export
conditional.tables <- function(x, variables, condition = NULL, condition.value = NULL, undef = NaN, ...) {
  #if (!order) margin <- sort.int(margin)
  n <- ntables(x)
  p <- tdim(x)
  if (max(variables, condition) > length(p)) stop("Not enough dimensions in table")
  if (!is.null(condition.value)) {
    dim_out <- c(p[variables],n,lengths(condition.value))
    rep_dim <- TRUE
  }
  else rep_dim <- FALSE
  if (is.list(condition.value)) condition.value = c(list(seq_len(n)), condition.value)
  out <- conditionTable(as.array(x), variables=variables+1, condition=c(1,condition+1), condition.value, undef=undef)
  
  ## in case conditionTable has dropped dimensions
  if (rep_dim) dim(out) <- dim_out
  out <- aperm.default(out, 
                       c(length(variables)+1, 
                         seq_along(variables), 
                         seq_along(condition)+length(variables)+1))
  
  
  attr(out, "tdim") <- dim(out)[-1]
  class(out) <- "tables"
  
  dim(out) = c(n, length(out)/n)
  
  out    
}

##' @method conditional2 tables
##' @export
conditional2.tables <- function(x, variables, condition = NULL, undef = NaN, ...) {
  n <- ntables(x)
  if (max(variables, condition) > length(tdim(x))) stop("Not enough dimensions in table")
  out <- conditionTable2(as.array(x), variables=variables+1, condition=c(1,condition+1), undef=undef)
  attributes(out)
  attr(out, "tdim") <- attr(x, "tdim")
  attr(out, "tdimnames") <- attr(x, "tdimnames")
  class(out) <- "tables"
  dim(out) = dim(x)
  
  out    
}

##' @method intervention default
##' @export
intervention.default <- function(x, variables, condition, ...) {
  interventionTable(x, variables, condition)
}

##' @method intervention tables
##' @export
intervention.tables <- function(x, variables, condition, ...) {
  tmp = conditional2(x, variables, condition, undef = .5)
  x = x/c(tmp)
  x
}


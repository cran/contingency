##' Calculates analytical bounds on the average causal effect from Balke and Pearl (1997)
##' 
##' @param x array containing a probability distribution, or object of class 
##' \code{tables} containing many probability distributions
##' @param instrument index of the instrument
##' @param treatment index of the treatment variable
##' @param outcome index of the outcome variable
##' 
##' @examples
##' set.seed(123)
##' p <- rprobdist(2,3)
##' acebounds(p, 1, 2, 3)
##' 
##' ps <- rprobMat(10,2,3)
##' acebounds(p, 1, 2, 3)
##' 
##' @export acebounds
acebounds <- function (x, instrument, treatment, outcome, ...) UseMethod("acebounds")

##' @export acebounds.default
acebounds.default <- function(x, instrument, treatment, outcome) {
  if (length(dim(x)) < 3) stop("Must be at least 3-dimensional array")
  tmp <- conditional(x, c(treatment, outcome), instrument)
  if (any(dim(tmp) != 2)) stop("Binary variables only")
  if (length(dim(tmp)) != 3) stop("Indices must be distinct")
  
  lwr = max(
    tmp[1,1,1]+tmp[2,2,2]-1,
    tmp[1,1,2]+tmp[2,2,2]-1,
    tmp[1,1,1]+tmp[2,2,1]-1,
    tmp[1,1,2]+tmp[2,2,1]-1,
    2*tmp[1,1,1]+tmp[2,2,1]+tmp[1,2,2]+tmp[2,2,2]-2,
    tmp[1,1,1]+2*tmp[2,2,1]+tmp[1,1,2]+tmp[2,1,2]-2,
    tmp[1,2,1]+tmp[2,2,1]+2*tmp[1,1,2]+tmp[2,2,2]-2,
    tmp[1,1,1]+tmp[2,1,1]+tmp[1,1,2]+2*tmp[2,2,2]-2)
  
  upp = min(1-tmp[1,2,1]-tmp[2,1,2],
            1-tmp[2,1,1]-tmp[1,2,2],
            1-tmp[2,1,1]-tmp[1,2,1],
            1-tmp[2,1,2]-tmp[1,2,2],
            2-2*tmp[2,1,1]-tmp[1,2,1]-tmp[1,2,2]-tmp[2,2,2],
            2-tmp[2,1,1]-2*tmp[1,2,1]-tmp[1,1,2]-tmp[2,1,2],
            2-tmp[2,2,1]-tmp[1,2,1]-tmp[1,2,2]-2*tmp[2,1,2],
            2-tmp[2,1,1]-tmp[1,1,1]-2*tmp[1,2,2]-tmp[2,1,2])
  c(lwr, upp)
}

##' Test instrumental inequalities
##' 
##' Check whether the instrumental inequalty is satisfied for 2x2x2 tables.
##' @param x object of class \code{tables}.
##' @param instrument index of the instrument
##' @param treatment index of the treatment variable
##' @param outcome index of the outcome variable
##' 
instIneq.tables <- function(x, instrument, treatment, outcome) {
  if (length(tdim(x)) < 3) stop("Must be at least 3-dimensional table")
  tmp <- conditional(x, c(treatment, outcome), instrument)
  if (any(tdim(tmp) != 2)) stop("Binary variables only")
  if (length(tdim(tmp)) != 3) stop("Indices must be distinct")
  
  pmax(tmp[,1,1,1]+tmp[,1,2,2],
       tmp[,2,1,1]+tmp[,2,2,2],
       tmp[,1,2,1]+tmp[,1,1,2],
       tmp[,2,2,1]+tmp[,2,1,2]) <= 1
}

acebounds.tables <- function(x, instrument, treatment, outcome) {
  if (length(tdim(x)) < 3) stop("Must be at least 3-dimensional table")
  tmp <- conditional(x, c(treatment, outcome), instrument)
  if (any(tdim(tmp) != 2)) stop("Binary variables only")
  if (length(tdim(tmp)) != 3) stop("Indices must be distinct")
  
  ## bounds from Cai et al (2008)
  lwr = pmax(
    tmp[,1,1,1]+tmp[,2,2,2]-1,
    tmp[,1,1,2]+tmp[,2,2,2]-1,
    tmp[,1,1,1]+tmp[,2,2,1]-1,
    tmp[,1,1,2]+tmp[,2,2,1]-1,
    2*tmp[,1,1,1]+tmp[,2,2,1]+tmp[,1,2,2]+tmp[,2,2,2]-2,
    tmp[,1,1,1]+2*tmp[,2,2,1]+tmp[,1,1,2]+tmp[,2,1,2]-2,
    tmp[,1,2,1]+tmp[,2,2,1]+2*tmp[,1,1,2]+tmp[,2,2,2]-2,
    tmp[,1,1,1]+tmp[,2,1,1]+tmp[,1,1,2]+2*tmp[,2,2,2]-2)
  
  upp = pmin(1-tmp[,1,2,1]-tmp[,2,1,2],
             1-tmp[,2,1,1]-tmp[,1,2,2],
             1-tmp[,2,1,1]-tmp[,1,2,1],
             1-tmp[,2,1,2]-tmp[,1,2,2],
             2-2*tmp[,2,1,1]-tmp[,1,2,1]-tmp[,1,2,2]-tmp[,2,2,2],
             2-tmp[,2,1,1]-2*tmp[,1,2,1]-tmp[,1,1,2]-tmp[,2,1,2],
             2-tmp[,2,2,1]-tmp[,1,2,1]-tmp[,1,2,2]-2*tmp[,2,1,2],
             2-tmp[,2,1,1]-tmp[,1,1,1]-2*tmp[,1,2,2]-tmp[,2,1,2])
  empty <- (lwr > upp)
  lwr[empty] = upp[empty] = NA

  cbind(lwr, upp)
}

##' Average Causal Effect bounds under weak constraints
##' 
##' The "message passing" implementation of constraint validation from WPP paper. 
##' Constraints are not the tightest possible, but don't require running a linear
##' programming solver.
##' 
##' @param p object of class \code{tables} containing 2x2x2 probability distributions
##' @param epsilons vector of relaxation parameters.  Defaults to IV model.
##' @param solutions logical: return full interval or just whether interval is non-empty?
##' 
##' @details If \code{solutions} is \code{TRUE} the function calculates intervals and 
##' returns them, giving \code{NA} if the distribution violates the constraints.
##' 
##' @return By default returns a matrix with two columns containing the lower and 
##' upper ACE bounds.  If \code{solutions = FALSE} returns logical vector with 
##' \code{TRUE} for those distributions which satisfy the constraints, \code{FALSE} 
##' otherwise.
##' 
##' @export aceBounds
##' @export aceBounds.default
##' @export aceBounds.tables
aceBounds <- function(p, instrument=1, treatment=2, outcome=3, epsilons=c(0,1,0,1,1,1), solutions = TRUE, ...) {
  UseMethod("aceBounds")
}

aceBounds.default <- function(p, instrument=1, treatment=2, outcome=3, epsilons=c(0,1,0,1,1,1), solutions = TRUE)
{
  p2 <- p
  dim(p2) <- c(1,prod(dim(p)))
  class(p2) <- "tables"
  tdim(p2) <- dim(p)
  out <- aceBounds.tables(p2, instrument, treatment, outcome, epsilons, solutions)[1,]
  
  return(out)
}

  
aceBounds.tables <- function(p, instrument=1, treatment=2, outcome=3, epsilons=c(0,1,0,1,1,1), solutions = TRUE)
{
  p_xy.w <- conditional(p, c(treatment, outcome), instrument)
  p_y.xw = conditional(p_xy.w, 2, c(1,3))
  p_x.w = conditional(p_xy.w, 1, 3)
  
  N <- ntables(p)
  pass <- rep(1, N)
  
  ## indices
  i_x.w = 1:4
  i_xp.w = c(2,1,4,3)
  i_x.wp = c(3,4,1,2)
  i_xp.wp = c(4,3,2,1)
  
  ## p_xy.w = matrix(c(P_YX.W0, P_YX.W1), nrow=N)
  
  eps_w      <- epsilons[1]
  eps_Y      <- epsilons[2]
  eps_X      <- epsilons[4]
  beta_lower <- epsilons[5]
  beta_upper <- epsilons[6]
  
  ## Calculate auxiliary variables
  UK_XY.W = pmin(p_xy.w / beta_lower, 1)
  LK_XY.W = p_xy.w / beta_upper
  
  UChi = beta_upper*p_x.w
  LChi = beta_lower*p_x.w
  
  L_X = pmax(p_x.w - eps_X, 0)
  U_X = pmin(p_x.w + eps_X, 1)
  
  U_Y = pmin(p_y.xw[,c(3,4,7,8), drop=FALSE] + eps_Y, 1)
  L_Y = pmax(p_y.xw[,c(3,4,7,8), drop=FALSE] - eps_Y, 0)
  
  U_bar <- rowMaxs(U_Y)
  L_bar <- rowMins(L_Y)
  
  ## Derive the box constraints for omega_xw first
  ## Theorem 1 bounds
  upper =                 UK_XY.W[,c(3,4,7,8), drop=FALSE] + U_Y * UChi[,i_xp.w, drop=FALSE]
  upper = pmin(upper,     UK_XY.W[,c(3,4,7,8), drop=FALSE] / L_X)
  upper = pmin(upper, 1 - LK_XY.W[,c(1,2,5,6), drop=FALSE] / U_X)
  
  lower =                 LK_XY.W[,c(3,4,7,8), drop=FALSE] + L_Y * LChi[,i_xp.w, drop=FALSE]
  lower = pmax(lower,     LK_XY.W[,c(3,4,7,8), drop=FALSE] / U_X)
  lower = pmax(lower, 1 - UK_XY.W[,c(1,2,5,6), drop=FALSE] / L_X)
  
  ## Theorem 2 bounds
  upper <- pmin(upper,     (UK_XY.W[,c(7,8,3,4), drop=FALSE] + eps_w * UChi[, i_x.wp, drop=FALSE]) / L_X[, i_x.wp, drop=FALSE])
  upper <- pmin(upper, 1 - (LK_XY.W[,c(5,6,1,2), drop=FALSE] - eps_w * UChi[, i_x.wp, drop=FALSE]) / U_X[, i_x.wp, drop=FALSE])
  lower <- pmax(lower,     (LK_XY.W[,c(7,8,3,4), drop=FALSE] - eps_w * UChi[, i_x.wp, drop=FALSE]) / U_X[, i_x.wp, drop=FALSE])
  lower <- pmax(lower, 1 - (UK_XY.W[,c(5,6,1,2), drop=FALSE] + eps_w * UChi[, i_x.wp, drop=FALSE]) / L_X[, i_x.wp, drop=FALSE])
  
  ## Bounds from Theorem 3
  upper <- pmin(upper,
                (UK_XY.W[,c(8,7,4,3), drop=FALSE] + UK_XY.W[,c(7,8,3,4), drop=FALSE] + UK_XY.W[,c(3,4,7,8), drop=FALSE] +
                   - LK_XY.W[,c(4,3,8,7), drop=FALSE] + UChi[,i_xp.w, drop=FALSE]*(U_bar + L_bar + 2 * eps_w) - L_bar))
  upper <- pmin(upper,
                (UK_XY.W[,c(4,3,8,7), drop=FALSE] + UK_XY.W[, c(7,8,3,4), drop=FALSE] + UK_XY.W[, c(3,4,7,8), drop=FALSE] - LK_XY.W[, c(8,7,4,3), drop=FALSE] +
                   + 2 * UChi[,i_xp.w, drop=FALSE] * eps_w + UChi[,c(4,3,2,1), drop=FALSE] * (U_bar + L_bar) - L_bar))
  lower <- pmax(lower,
                (- UK_XY.W[,c(8,7,4,3), drop=FALSE] + LK_XY.W[,c(7,8,3,4), drop=FALSE] + LK_XY.W[,c(3,4,7,8), drop=FALSE] +
                   LK_XY.W[,c(4,3,8,7), drop=FALSE] + - 2 * UChi[,i_xp.w, drop=FALSE] * eps_w + LChi[,c(4,3,2,1), drop=FALSE] * (U_bar + L_bar) - U_bar))
  lower <-  pmax(lower,
                 (-UK_XY.W[,c(4,3,8,7), drop=FALSE] + LK_XY.W[,c(7,8,3,4), drop=FALSE] + LK_XY.W[,c(3,4,7,8), drop=FALSE] +
                    + LK_XY.W[,c(8,7,4,3), drop=FALSE] - UChi[,i_xp.w, drop=FALSE] * 2 * eps_w + LChi[,i_xp.w, drop=FALSE]*(U_bar +L_bar) - U_bar))
  
  upper[is.nan(upper)] = 1
  upper = pmin(upper, 1)
  lower[is.nan(lower)] = 0
  lower = pmax(lower, 0)
  
  pass = pass*rowMins(upper >= lower)
  
  ## bounds for omega_xw
  omega_upper = upper
  omega_lower = lower
  
  ## bounds for differences omega_xw - omega_xw'
  diff_upper <- matrix(eps_w, nrow=N, ncol = 4)
  diff_lower <- matrix(-eps_w, nrow=N, ncol = 4)
  
  ## Now, iterate over linear constraints
  for (iter in 1:4) {    
    upper = omega_upper
    lower = omega_lower
    
    ## #############################################
    ## Iteration over the linear constraints of Theorem 2
    upper = pmin(upper,
                 omega_upper[,i_x.wp, drop=FALSE]*U_X[,i_xp.w, drop=FALSE] +
                   UK_XY.W[,c(3,4,7,8), drop=FALSE] + eps_w*UChi[,i_xp.w, drop=FALSE])
    upper = pmin(upper,
                 (omega_upper[,i_x.wp, drop=FALSE] - 1)*L_X[,i_xp.w, drop=FALSE] +
                   1 - LK_XY.W[,c(1,2,5,6), drop=FALSE]  + eps_w * UChi[,i_xp.w, drop=FALSE])
    upper = pmin(upper, omega_upper[,i_x.wp, drop=FALSE] + eps_w)
    
    lower = pmax(lower,
                 omega_lower[,i_x.wp, drop=FALSE]*L_X[,i_xp.w, drop=FALSE] +
                   LK_XY.W[,c(3,4,7,8), drop=FALSE] - eps_w * UChi[,i_xp.w, drop=FALSE])
    lower = pmax(lower,
                 (omega_lower[,i_x.wp, drop=FALSE] - 1)*U_X[,i_xp.w, drop=FALSE] +
                   1 - UK_XY.W[,c(1,2,5,6), drop=FALSE] - eps_w * UChi[,i_xp.w, drop=FALSE])
    lower = pmax(lower, omega_lower[,i_x.wp, drop=FALSE] - eps_w)
    
    omega_upper = upper
    omega_lower = lower
    
    pass = pass*rowMins(upper >= lower)
    
    ## #############################################
    ## Iteration over the linear constraints of Theorem 2 to bound omega_xw - omega_xw'
    ## equation (9)
    upper = pmin(diff_upper,
                 omega_upper[,i_x.wp, drop=FALSE]*(U_X[,i_xp.w, drop=FALSE] - 1) +
                   UK_XY.W[,c(3,4,7,8), drop=FALSE] + eps_w*UChi[,i_xp.w, drop=FALSE])
    upper = pmin(upper,
                 (omega_upper[,i_x.wp, drop=FALSE] - 1)*(L_X[,i_xp.w, drop=FALSE] - 1)  +
                   - LK_XY.W[,c(1,2,5,6), drop=FALSE]  + eps_w * UChi[,i_xp.w, drop=FALSE])
    upper = pmin(upper, omega_upper - omega_lower[,i_x.wp, drop=FALSE])
    
    lower = pmax(diff_lower,
                 omega_lower[,i_x.wp, drop=FALSE]*(L_X[,i_xp.w, drop=FALSE]-1) +
                   LK_XY.W[,c(3,4,7,8), drop=FALSE] - eps_w * UChi[,i_xp.w, drop=FALSE])
    lower = pmax(lower,
                 (omega_lower[,i_x.wp, drop=FALSE] - 1)*(U_X[,i_xp.w, drop=FALSE]-1) +
                   - UK_XY.W[,c(1,2,5,6), drop=FALSE] - eps_w * UChi[,i_xp.w, drop=FALSE])
    lower = pmax(lower, omega_lower - omega_upper[,i_x.wp, drop=FALSE])
    
    ## ## in the old code we did the following, but this doesn't look right
    ## upper = pmax(upper, 0)
    ## lower = pmax(lower, 0)
    
    diff_upper = upper
    diff_lower = lower
    
    ## #############################################
    ## Iteration over the linear constraints of Theorem 3 to bound omega_xw
    upper = pmin(omega_upper,
                 diff_upper[,i_xp.w, drop=FALSE] - LK_XY.W[,c(4,3,8,7), drop=FALSE] + UK_XY.W[,c(3,4,7,8), drop=FALSE] +
                   + UK_XY.W[,c(8,7,4,3), drop=FALSE] + UK_XY.W[,c(7,8,3,4), drop=FALSE] +
                   - LChi[,i_x.w, drop=FALSE]*(U_bar + L_bar) + 2*eps_w + UChi[,i_x.wp, drop=FALSE] + U_bar)
    upper = pmin(upper,
                 diff_upper[,i_xp.wp, drop=FALSE] - LK_XY.W[,c(8,7,4,3), drop=FALSE] + UK_XY.W[,c(7,8,3,4), drop=FALSE] +
                   + UK_XY.W[,c(4,3,8,7), drop=FALSE] + UK_XY.W[,c(3,4,7,8), drop=FALSE] +
                   + 2*eps_w*UChi[,i_x.wp, drop=FALSE] - LChi[,i_x.wp, drop=FALSE]*(U_bar + L_bar) + U_bar)
    
    lower = pmax(omega_lower,
                 - diff_lower[,i_xp.wp, drop=FALSE] - UK_XY.W[,c(8,7,4,3), drop=FALSE] + LK_XY.W[,c(3,4,7,8), drop=FALSE] +
                   + LK_XY.W[,c(4,3,8,7), drop=FALSE] + LK_XY.W[,c(7,8,3,4), drop=FALSE] +
                   - UChi[,i_x.wp, drop=FALSE]*(U_bar + L_bar + 2*eps_w) + L_bar)
    lower = pmax(lower,
                 - diff_lower[,i_xp.w, drop=FALSE] - UK_XY.W[,c(4,3,8,7), drop=FALSE] + LK_XY.W[,c(7,8,3,4), drop=FALSE] +
                   + LK_XY.W[,c(8,7,4,3), drop=FALSE] + LK_XY.W[,c(3,4,7,8), drop=FALSE] +
                   - 2*eps_w*UChi[,i_x.wp, drop=FALSE] - UChi[,i_x.w, drop=FALSE]*(U_bar + L_bar) + L_bar)
    
    omega_upper = upper
    omega_lower = lower
    
    pass = pass*rowMins(omega_upper >= omega_lower)
  }
  
  if (solutions)  {
    intervals <- matrix(0, nrow = N, ncol = 2)
    
    alpha_upper = beta_upper * pmin(omega_upper, 1)
    alpha_lower = beta_lower * omega_lower
    
    p_w <- conditional(p, instrument, c())
    intervals[, 1] <- (alpha_lower[, 4] - alpha_upper[, 3]) * p_w[,2] + (alpha_lower[, 2] - alpha_upper[, 1]) * p_w[,1]
    intervals[, 2] <- (alpha_upper[, 4] - alpha_lower[, 3]) * p_w[,2] + (alpha_upper[, 2] - alpha_lower[, 1]) * p_w[,1]
#    intervals[pass == 0, 1] <- NA
#    intervals[pass == 0, 2] <- NA
    return(intervals)
  }
  
  return(pass)
}

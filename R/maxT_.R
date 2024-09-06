
#' @title Westfall & Young's max\eqn{T} Algorithm
#' 
#' @description
#' Westfall & Young's max\eqn{T} algorithm, as described in Box 2, page 82 of \doi{10.1214/ss/1056397487}
#' 
#' @param t. \link[base]{numeric} \link[base]{vector} of 
#' \eqn{t}-statistic \eqn{(t_1,\cdots,t_m)^t} 
#' for each hypothesis \eqn{H_j, j = 1,\cdots,m}
#' 
#' @param T. \link[base]{numeric} \link[base]{matrix} of dimension \eqn{(m,B)},
#' permuted \eqn{t}-statistic \eqn{(t_{1,b},\cdots,t_{m,b})^t}
#' for permutation \eqn{b=1,\cdots,B}
#' 
#' @param two.sided \link[base]{logical} scalar,
#' whether to take \link[base]{abs}olute of `t.` and `T.`
#' before other calculations.
#' Default `TRUE` as in \doi{10.1214/ss/1056397487},
#' while Moodie's algorithm (\url{https://rundfr.fredhutch.org}) does not have this step.
#' 
#' @references
#' Peter H. Westfall, S. Stanley Young (1993). *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. \url{https://www.wiley.com/en-us/Resampling-Based+Multiple+Testing%3A+Examples+and+Methods+for+p-Value+Adjustment-p-9780471557616}
#' 
#' S. Dudoit, J. P. Shaffer, J. C. Boldrick (2003). *Multiple Hypothesis Testing in Microarray Experiments*, Statistical Science 18(1), 71-103 
#' \doi{10.1214/ss/1056397487}
#' 
#' 
#' @export
maxT_ <- function(t., T., two.sided = TRUE) {
  
  m <- length(t.)
  if (!is.matrix(T.)) stop('`T.` must be matrix')
  dm <- dim(T.)
  if (dm[1L] != m) stop('nrow(T.) need to be length(t.)')
  if (anyNA(t.) || anyNA(T.)) stop('do not allow NA in `t.` or `T.`')
  
  # I am not sure if parameter `two.sided` is appropriate
  # Moodie's code does *not* have ?base::abs
  # \doi{10.1214/ss/1056397487} has ?base::abs
  if (two.sided) {
    t. <- abs(t.)
    T. <- abs(T.)
  }
  
  r <- order(t., decreasing = TRUE)
  # \eqn{(r_1, \cdots, r_m)^t}
  
  # ordered `t.`
  tr <- t.[r] 
  # \eqn{|t_{r_1}| \geq \cdots \geq |t_{r_m}|}
  
  # Step 3: successive maxima
  # `U`: matrix of dimension \eqn{(m,B)}
  U <- if (m == 1L) T. else {
    # order rows of `T.` by `r`
    Tr <- T.[r, , drop = FALSE] # \eqn{|t_{r_1,b}|, \cdots, |t_{r_m,b}|}
    apply(Tr, MARGIN = 2L, FUN = function(x) {
      # 'successive maxima'
      rev.default(cummax(rev.default(x)))
      # \eqn{u_{m,b}} and recursive \eqn{u_{j,b}}
    })
  }
  
  # permutation adjusted p-values
  # \eqn{\tilde{p}^*_{r_j}}
  p0 <- .rowMeans(U >= tr, m = m, n = dm[2L], na.rm = TRUE)
  
  # monotonicity constraints enforced
  p <- cummax(p0)
  
  # restore the order of `t.`
  return(p[order(r)])
  
}


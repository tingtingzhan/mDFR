

#' @title Westfall & Young's \linkS4class{maxT} Algorithm
#' 
#' @description
#' Westfall & Young's \linkS4class{maxT} algorithm, as described in Box 2, page 82 of \doi{10.1214/ss/1056397487}.
#' 
#' @note
#' The algorithm is implemented in an unexported
#' \link[methods]{initialize} method, which could be revealed by 
#' `getMethod(initialize, signature = 'maxT')` to curious eyes.
#' 
#' 
#' @slot t. \link[base]{length}-\eqn{m} \link[base]{double} \link[base]{vector},
#' test statistics \eqn{t_1,\cdots,t_m} 
#' for each hypothesis \eqn{H_j}, \eqn{j = 1,\cdots,m},
#' in the original data
#' 
#' @slot T. \link[base]{double} \link[base]{matrix} of \link[base]{dim}ension \eqn{(m,B)},
#' test statistics \eqn{t_{j,b}}
#' for permutation \eqn{b=1,\cdots,B}
#' 
#' @slot tr \link[base]{length}-\eqn{m} \link[base]{double} \link[base]{vector}, 
#' ordered test statistics 
#' \eqn{t_{r_1}\geq\cdots\geq t_{r_m}} for one-sided test, or 
#' \eqn{|t_{r_1}|\geq\cdots\geq|t_{r_m}|} for two-sided test
#' 
#' @slot U \link[base]{double} \link[base]{matrix} of \link[base]{dim}ension \eqn{(m,B)}, 
#' the successive maxima \eqn{u_{j,b}}
#' 
#' @slot p_perm \link[base]{length}-\eqn{m} \link[base]{double} \link[base]{vector}, 
#' permutation adjusted \eqn{p}-values \eqn{p_{r_j}}
#' 
#' @slot p_mono \link[base]{length}-\eqn{m} \link[base]{double} \link[base]{vector}, 
#' permutation adjusted \eqn{p}-values under monotonicity constraints \eqn{\tilde{p}_{r_j}}
#' 
#' @slot p. \link[base]{length}-\eqn{m} \link[base]{double} \link[base]{vector}, 
#' \eqn{\tilde{p}_{r_j}} restored in the order of original hypotheses \eqn{1,\cdots,m}
#' 
#' @slot alternative \link[base]{character} scalar, either 
#' `'two.sided'` (default value), `'greater'` or `'less'`
#' 
#' @slot design (optional) \link[base]{data.frame}, 
#' the study design
#' 
#' @slot label (optional) \link[base]{character} scalar, 
#' the study name
#' 
#' @references
#' Peter H. Westfall, S. Stanley Young (1993). *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. \url{https://www.wiley.com/en-us/Resampling-Based+Multiple+Testing%3A+Examples+and+Methods+for+p-Value+Adjustment-p-9780471557616}
#' 
#' S. Dudoit, J. P. Shaffer, J. C. Boldrick (2003). *Multiple Hypothesis Testing in Microarray Experiments*, 
#' \doi{10.1214/ss/1056397487}
#' 
#' \url{https://tingtingzhan-maxt.netlify.app/appendix/westfall_yang.html}
#' 
# @name maxT
# @aliases maxT-class
#' @export
setClass(Class = 'maxT', slots = c(
  t. = 'numeric', T. = 'matrix',
  tr = 'numeric', U = 'matrix',
  p_perm = 'numeric', p_mono = 'numeric', 
  p. = 'numeric',
  alternative = 'character',
  design = 'data.frame',
  label = 'character'
), prototype = prototype(
  alternative = 'two.sided'
))


# @param .Object the \linkS4class{maxT} object to be \link[methods]{initialize}d
#' 
# @param ... additional parameters, currently not in use
#' 

setMethod(f = initialize, signature = 'maxT', definition = \(.Object, ...) {
  
  x <- callNextMethod(.Object, ...)
  
  t. <- x@t. 
  T. <- x@T.

  m <- length(t.)
  if (!is.matrix(T.)) stop('`T.` must be matrix')
  dm <- dim(T.)
  if (dm[1L] != m) stop('nrow(T.) need to be length(t.)')
  
  if (anyNA(T.)) stop('Does not allow missingness in `T.`')
  if (anyNA(t.)) stop('Does not allow missingness in `t.`')
  
  switch(x@alternative, 'two.sided' = {
    t. <- abs(t.)
    T. <- abs(T.)
  }, 'greater' = {
    # do nothing
  }, 'less' = {
    t. <- - t.
    T. <- - T.
  })
  
  o <- order(t., decreasing = TRUE)
  # \eqn{(r_1, \cdots, r_m)^t}
  
  # ordered `t.`
  tr <- t.[o] 
  
  # successive maxima
  # `U`: matrix of dimension \eqn{(m,B)}
  U <- if (m == 1L) T. else {
    # order rows of `T.` by `o`
    Tr <- T.[o, , drop = FALSE] # \eqn{|t_{r_1,b}|, \cdots, |t_{r_m,b}|}
    apply(Tr, MARGIN = 2L, FUN = \(x) {
      # 'successive maxima'
      x |>
        rev.default() |>
        cummax() |>
        rev.default()
      # \eqn{u_{m,b}} and recursive \eqn{u_{j,b}}
    })
  }
  
  # permutation adjusted p-values
  # \eqn{p_{r_j}}
  p_perm <- .rowMeans(U >= tr, m = m, n = dm[2L], na.rm = TRUE)
  
  # monotonicity constraints enforced
  # \eqn{\tilde{p}_{r_j}}
  p_mono <- cummax(p_perm)
  
  x@tr <- tr 
  x@U <- U
  x@p_perm <- p_perm
  x@p_mono <- p_mono
  x@p. <- p_mono[order(o)]
  return(x)
  
})







#' @method as.data.frame maxT
#' @export
as.data.frame.maxT <- function(x, ..., check.names = FALSE) {
  
  tmp <- list(
    
    if (length(x@design)) x@design, #  else NULL
    
    't' = if (x@alternative == 'greater') {
      x@t. |> 
        round(digits = 3L)
    }, # else NULL
    
    'negative.t' = if (x@alternative == 'less') {
      (- x@t.) |> 
        round(digits = 3L)
    }, # else NULL
    
    'abs(t)' = if (x@alternative == 'two.sided') {
      abs(x@t.) |> 
        round(digits = 3L)
    }, # else NULL
    
    p.adj = x@p. |> 
      label_pvalue_sym()()
    
  )
  
  tmp[lengths(tmp) > 0L] |>
    as.data.frame.list(check.names = FALSE, ...)
  
}












setMethod(f = show, signature = 'maxT', definition = \(object) {
  object |>
    as.data.frame.maxT() |>
    print()
})







#' @importFrom ggplot2 autoplot ggplot
#' @export
autoplot.maxT <- function(object, ...) {
  suppressWarnings(expr = {
    ggplot() + autolayer.maxT(object, ...)
  })
  # Vectorized input to `element_text()` is not officially supported.
  # Results may be unexpected or may change in future versions of ggplot2. 
}



# @param conf.level \link[base]{double} scalar, 
# confidence level, or \eqn{1-\alpha}. Default .95 (or \eqn{\alpha=.05})
# @details
# The function [autoplot.maxT()] plots 
# \itemize{
# \item the successive maxima 
# \eqn{u_{jb}}, \eqn{j=1,\cdots,m}, \eqn{b=1,\cdots,B}, and 
# \item the decreasing-ordered statistics \eqn{|t_{r_1}|\geq|t_{r_2}|\geq\cdots\geq|t_{r_m}|} for two-sided test,
# or \eqn{t_{r_1}\geq t_{r_2}\geq\cdots\geq t_{r_m}} for one-sided test
# }
# Printed on opposing axis are
# \itemize{
# \item values of the decreasing-ordered statistics
# \item permutation adjusted \eqn{p}-values \eqn{\tilde{p}_{r_j}}, as well as under monotonicity constraints \eqn{\tilde{p}^*_{r_j}}
# }
# Tests with \eqn{\tilde{p}_{r_j}\leq\alpha} is considered significant
# and colored pink (hex color `#F8766D`), otherwise non-significant and colored blue (hex color `#00BFC4`)
#' @importFrom ggplot2 autolayer aes geom_jitter geom_point geom_line sec_axis scale_x_continuous scale_y_continuous labs theme element_text
#' @importFrom scales pal_hue
#' @export
autolayer.maxT <- function(object, conf.level = .95, ...) {

  p_perm <- object@p_perm 
  p_mono <- object@p_mono
  tr <- object@tr
  U <- object@U
  
  dm <- dim(U)
  
  col0 <- pal_hue()(n = 2L)
  col <- ifelse(p_mono <= (1 - conf.level), yes = col0[1L], no = col0[2L])
  
  tseq <- seq_along(tr)
  
  mp_point <- aes(y = tseq, x = tr)
  
  list(
    geom_jitter(
      mapping = aes(y = rep(tseq, each = dm[2L]), x = c(t.default(U))), 
      color = rep(col, each = dm[2L]), 
      width = 0, height = .25, # no need to jitter on width!
      size = 1e-5, alpha = .1, show.legend = FALSE
    ),
    
    geom_point(mapping = mp_point, color = col, size = 1.5, show.legend = FALSE),
    
    geom_line(mapping = mp_point, color = rev.default(col), linewidth = .5, show.legend = FALSE),
    
    scale_y_continuous(
      name = 'Permutation Adjusted p-values \u27a4 Monotonicity Constraints',
      breaks = tseq, 
      minor_breaks = NULL, 
      labels = sprintf(fmt = '%.3f \u27a4 %.3f', p_perm, p_mono),
      sec.axis = sec_axis(
        transform = ~.,
        name = 'Test-Statistic',
        breaks = tseq,
        labels = tr |> 
          sprintf(fmt = '%.2f') # |> 
          # format(justify = 'right') # no need!
      )
    ),
    
    theme(axis.text.y = element_text(colour = col)),
    
    labs(
      x = 'Successive Maxima of Permuted Test-Statistic'
    )
  )
  
}




#' @export
'[.maxT' <- function(x, i) {
  new(
    Class = 'maxT', 
    t. = x@t.[i], 
    T. = x@T.[i, , drop = FALSE], 
    two.sided = x@two.sided
  )
}



setOldClass('gtable')

#' @title Westfall & Young's \linkS4class{maxT} Algorithm
#' 
#' @description
#' Westfall & Young's \linkS4class{maxT} algorithm, as described in Box 2, page 82 of \doi{10.1214/ss/1056397487}.
#' 
# @param .Object the \linkS4class{maxT} object to be \link[methods]{initialize}d
#' 
# @param ... additional parameters, currently not in use
#' 
#' @details
#' In the original data, obtain the test statistics \eqn{t_1,\cdots,t_m} 
#' for each hypothesis \eqn{H_j}, \eqn{j=1,\cdots,m},
#' as well as the *decreasing order statistics* \eqn{\mathbf{r} = (r_1,\cdots,r_m)^t}, such that 
#' \deqn{
#' \begin{cases}
#' |t_{r_1}|\geq|t_{r_2}|\geq\cdots\geq|t_{r_m}| & \text{for two-sided test}\\
#' t_{r_1}\geq t_{r_2}\geq\cdots\geq t_{r_m} & \text{for one-sided test}
#' \end{cases}
#' } 
#' 
#' In each permuted copy \eqn{b}, \eqn{b=1,\cdots,B},
#' obtain the test statistics \eqn{t_{1,b},\cdots,t_{m,b}}
#' and their **successive maxima** for two-sided test
#' \deqn{
#' \begin{cases}
#' u_{m,b} & = |t_{r_m,b}|\\
#' u_{j,b} & = \text{max}\left\{u_{j+1,b}, |t_{r_j,b}|\right\}, \quad \text{for}\ j = m−1,\cdots,1\\
#' \end{cases}
#' }
#' or for one-sided test
#' \deqn{
#' \begin{cases}
#' u_{m,b} & = t_{r_m,b}\\
#' u_{j,b} & = \text{max}\left\{u_{j+1,b}, t_{r_j,b}\right\}, \quad \text{for}\ j = m−1,\cdots,1\\
#' \end{cases}
#' }
#' 
#' The permutation adjusted \eqn{p}-values are
#' \deqn{
#' \tilde{p}_{r_j} = 
#' \begin{cases}
#' \dfrac{1}{B}\displaystyle\sum_{b=1}^B I\left(u_{j,b}\geq|t_{r_j}|\right) & \text{for two-sided test}\\
#' \dfrac{1}{B}\displaystyle\sum_{b=1}^B I\left(u_{j,b}\geq t_{r_j}\right) & \text{for one-sided test}
#' \end{cases}
#' }
#' 
#' To enforce the monotonicity constraints,
#' \deqn{
#' \begin{cases}
#' \tilde{p}^*_{r_1} & = \tilde{p}_{r_1}\\
#' \tilde{p}^*_{r_j} & = \text{max}\left\{\tilde{p}_{r_j}, \tilde{p}^*_{r_j-1}\right\}, \quad \text{for}\ j = 2,\cdots,m\\
#' \end{cases}
#' }
#' 
#' @note
#' The altorithm described in **Details** is implemented in an unexported
#' \link[methods]{initialize} method, which could be revealed by 
#' `getMethod(initialize, signature = 'maxT')` to curious eyes.
#' 
#' 
#' @slot t. \link[base]{double} \link[base]{vector},
#' test statistics \eqn{t_1,\cdots,t_m} 
#' for each hypothesis \eqn{H_j}, \eqn{j = 1,\cdots,m},
#' in the original data
#' 
#' @slot T. \link[base]{double} \link[base]{matrix} of dimension \eqn{(m,B)},
#' test statistics \eqn{t_{j,b}}, \eqn{j = 1,\cdots,m},
#' for permutation \eqn{b=1,\cdots,B}
#' 
#' @slot tr \link[base]{double} \link[base]{vector}, ordered test statistics \eqn{t_{r_1}\geq t_{r_2}\geq\cdots\geq t_{r_m}} for one-sided test, or \eqn{|t_{r_1}|\geq|t_{r_2}|\geq\cdots\geq|t_{r_m}|} for two-sided test
#' @slot U \link[base]{double} \link[base]{matrix} of dimension \eqn{(m,B)}, successive maxima \eqn{u_{j,b}}, \eqn{j=1,\cdots,m}, \eqn{b=1,\cdots,B}
#' @slot p_perm \link[base]{double} \link[base]{vector}, permutation adjusted \eqn{p}-values \eqn{\tilde{p}_{r_j}}
#' @slot p_mono \link[base]{double} \link[base]{vector}, permutation adjusted \eqn{p}-values under monotonicity constraints \eqn{\tilde{p}^*_{r_j}}
#' @slot p. \link[base]{double} \link[base]{vector}, \eqn{\tilde{p}^*_{r_j}} restored in the order of original test statistics \eqn{t_1,\cdots,t_m}
#' 
#' @slot two.sided two.sided \link[base]{logical} scalar,
#' whether to perform a two sided test.
#' Default `TRUE` as in \doi{10.1214/ss/1056397487},
#' while Moodie's (\url{https://rundfr.fredhutch.org}) use `FALSE`
#' 
#' @slot design (optional) \link[base]{data.frame}, study design
#' @slot name (optional) \link[base]{character} scalar, study name
#' @slot gtable returned object of function \link[ggplot2]{ggplot_gtable}
#' 
#' @references
#' Peter H. Westfall, S. Stanley Young (1993). *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. \url{https://www.wiley.com/en-us/Resampling-Based+Multiple+Testing%3A+Examples+and+Methods+for+p-Value+Adjustment-p-9780471557616}
#' 
#' S. Dudoit, J. P. Shaffer, J. C. Boldrick (2003). *Multiple Hypothesis Testing in Microarray Experiments*, 
#' \doi{10.1214/ss/1056397487}
#' 
# @name maxT
# @aliases maxT-class
#' @importFrom ggplot2 ggplot_build ggplot_gtable
#' @importFrom methods setClass
#' @export
setClass(Class = 'maxT', slots = c(
  t. = 'numeric', T. = 'matrix',
  tr = 'numeric', U = 'matrix',
  p_perm = 'numeric', p_mono = 'numeric', 
  p. = 'numeric',
  two.sided = 'logical',
  design = 'data.frame',
  name = 'character',
  gtable = 'gtable'
), prototype = prototype(
  two.sided = TRUE,
  gtable = ggplot() |>
    ggplot_build() |>
    ggplot_gtable()
))


#' @importFrom ggplot2 ggplot_build ggplot_gtable
#' @importFrom methods setMethod callNextMethod initialize
setMethod(f = initialize, signature = 'maxT', definition = function(.Object, ...) {
  
  x <- callNextMethod(.Object, ...)
  
  t. <- x@t. # t.orig <- 
  T. <- x@T. # T.orig <- 
  two.sided <- x@two.sided
  
  m <- length(t.)
  if (!is.matrix(T.)) stop('`T.` must be matrix')
  dm <- dim(T.)
  if (dm[1L] != m) stop('nrow(T.) need to be length(t.)')
  
  if (anyNA(T.)) stop('Does not allow missingness in `T.`')
  if (anyNA(t.)) stop('Does not allow missingness in `t.`')
  
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
    apply(Tr, MARGIN = 2L, FUN = \(x) {
      # 'successive maxima'
      rev.default(cummax(rev.default(x)))
      # \eqn{u_{m,b}} and recursive \eqn{u_{j,b}}
    })
  }
  
  # permutation adjusted p-values
  # \eqn{\tilde{p}_{r_j}}
  p_perm <- .rowMeans(U >= tr, m = m, n = dm[2L], na.rm = TRUE)
  
  # monotonicity constraints enforced
  # \eqn{\tilde{p}^*_{r_j}}
  p_mono <- cummax(p_perm)
  
  # restore the order of `t.`
  if (FALSE) {
    set.seed(13134); x = rnorm(10L)
    o = order(x)
    s = sort(x)
    identical(x[order(x)], s)
    # Q: if I have `s` and `o`, how to restore `x`?
    identical(s[order(o)], x) # solution!
  } # R intro 101
  
  # `U` can be huge
  # rendering `gg` can be slow
  # read carefully!!!
  # https://stackoverflow.com/questions/73470828/ggplot2-is-slow-where-is-the-bottleneck
  # save slot `@gtable` inside \linkS4class{maxT}, to save time in printing
  
  gtb <- autoplot_maxT_(p_perm = p_perm, p_mono = p_mono, tr = tr, U = U, two.sided = two.sided) |> # always fast
    ggplot_build() |> # 'ggplot_built'
    ggplot_gtable() # c('gtable', 'gTree', 'grob', 'gDesc')
  # ggplot_build() and ggplot_gtable() not too slow
  # use parallel::mclapply() in production code!
  
  x@tr <- tr 
  x@U <- U
  x@p_perm <- p_perm
  x@p_mono <- p_mono
  x@p. <- p_mono[order(r)]
  x@gtable <- gtb
  return(x)
})







#' @title Convert \linkS4class{maxT} to a \link[base]{data.frame}
#' 
#' @param x a \linkS4class{maxT} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' Function [as.data.frame.maxT] returns a \link[base]{data.frame}.
#' 
#' @method as.data.frame maxT
#' @export as.data.frame.maxT
#' @export
as.data.frame.maxT <- function(x, ...) {
  tmp <- list(
    if (length(x@design)) x@design, #  else NULL
    tstat = x@t.,
    'abs(tstat)' = if (x@two.sided) abs(x@t.), # else NULL
    adjp = x@p.
  )
  as.data.frame.list(tmp[lengths(tmp) > 0L], check.names = FALSE)
}



#' @importFrom methods setMethod show 
setMethod(f = show, signature = 'maxT', definition = function(object) {
  object |>
    reactable_maxT() |>
    print() # ?htmlwidgets:::print.htmlwidget
})




#' @title Visualize Westfall & Young's \linkS4class{maxT} Algorithm
#' 
#' @description
#' To visualize Westfall & Young's \linkS4class{maxT} algorithm using package \CRANpkg{ggplot2}.
#' 
#' @param object a \linkS4class{maxT} object
#' 
#' @param p_perm,p_mono,tr,U,two.sided description
#' 
#' @param conf.level \link[base]{double} scalar, 
#' confidence level, or \eqn{1-\alpha}. Default .95 (or \eqn{\alpha=.05})
#' 
#' @param vertical \link[base]{logical} scalar, whether the
#' labels of ordered test statistics and permutation adjusted \eqn{p}-values should be on vertical axes,
#' default `TRUE`
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' Function [autoplot.maxT] plots 
#' \itemize{
#' \item the successive maxima 
#' \eqn{u_{jb}}, \eqn{j=1,\cdots,m}, \eqn{b=1,\cdots,B}, and 
#' \item the decreasing-ordered statistics \eqn{|t_{r_1}|\geq|t_{r_2}|\geq\cdots\geq|t_{r_m}|} for two-sided test,
#' or \eqn{t_{r_1}\geq t_{r_2}\geq\cdots\geq t_{r_m}} for one-sided test
#' }
#' Printed on opposing axis are
#' \itemize{
#' \item values of the decreasing-ordered statistics
#' \item permutation adjusted \eqn{p}-values \eqn{\tilde{p}_{r_j}}, as well as under monotonicity constraints \eqn{\tilde{p}^*_{r_j}}
#' }
#' Tests with \eqn{\tilde{p}^*_{r_j}\leq\alpha} is considered significant
#' and colored pink (hex color `#F8766D`), otherwise non-significant and colored blue (hex color `#00BFC4`)
#' 
#' See full details of these notations in \linkS4class{maxT}.
#' 
#' @returns
#' Function [autoplot.maxT] returns a \link[ggplot2]{ggplot} object.
#' 
#' @importFrom ggplot2 autoplot ggplot aes geom_jitter geom_point geom_line sec_axis scale_x_continuous scale_y_continuous labs theme
#' @importFrom ggtext element_markdown
#' @importFrom scales pal_hue
#' @name autoplot.maxT
#' @export autoplot.maxT
#' @export
autoplot.maxT <- function(object, ...) {
  autoplot_maxT_(p_perm = object@p_perm, p_mono = object@p_mono, tr = object@tr, U = object@U, two.sided = object@two.sided, ...)
}

#' @rdname autoplot.maxT
#' @export
autoplot_maxT_ <- function(
    p_perm, p_mono, tr, U, two.sided,
    conf.level = .95, 
    vertical = TRUE, 
    ...
) {
  
  dm <- dim(U)
  
  col0 <- pal_hue()(n = 2L)
  col <- ifelse(p_mono <= (1-conf.level), yes = col0[1L], no = col0[2L])
  
  tseq <- seq_along(tr)
  
  if (vertical) {
    mp_jitter <- aes(y = rep(tseq, each = dm[2L]), x = c(t.default(U)))
    mp_point <- aes(y = tseq, x = tr)
    scale_t_p <- scale_y_continuous 
  } else {
    mp_jitter <- aes(x = rep(tseq, each = dm[2L]), y = c(t.default(U)))
    mp_point <- aes(x = tseq, y = tr)
    scale_t_p <- scale_x_continuous
  }
  
  fig <- ggplot() + 
    geom_jitter(mapping = mp_jitter, color = rep(col, each = dm[2L]), width = .25, height = .25, size = 1e-5, alpha = .1, show.legend = FALSE) + 
    geom_point(mapping = mp_point, color = col, size = 1.5, show.legend = FALSE) + 
    geom_line(mapping = mp_point, color = if (vertical) rev.default(col) else col, linewidth = .5, show.legend = FALSE) +
    scale_t_p(
      name = 'Permutation Adjusted p-values (w. Monotonicity Constraints)',
      breaks = tseq, 
      minor_breaks = NULL, 
      labels = sprintf(fmt = '%.2f%% (%.2f%%)', 1e2*p_perm, 1e2*p_mono), # '\u27a4' is a beautiful arrow
      sec.axis = sec_axis(
        transform = ~.,
        name = sprintf(fmt = 'Original Test-Statistic, %s Sided', ifelse(two.sided, 'Two', 'One')),
        breaks = tseq,
        labels = sprintf(fmt = '%.2f', tr)
      )
    )
  
  U_lab <- sprintf(fmt = 'Successive Maxima of Permuted Test-Statistic, %s Sided', ifelse(two.sided, 'Two', 'One'))
  
  # https://github.com/kassambara/survminer/pull/503
  # 'the function ?ggtext::element_markdown is now used in place of ?ggplot2::element_text to handle vectorized colors'
  
  if (vertical) {
    fig + labs(x = U_lab) + theme(
      axis.text.y.left = element_markdown(colour = col),
      axis.text.y.right = element_markdown(colour = col)
    )
  } else {
    fig + labs(y = U_lab) + theme(
      axis.text.x.top = element_markdown(angle = 90, vjust = 0.5, color = col),
      axis.text.x.bottom = element_markdown(angle = 90, vjust = 0.5, color = col)
    )
  }
    
}



#' @title Subset of Hypothesis in Westfall & Young's \linkS4class{maxT} Algorithm
#' 
#' @param x a \linkS4class{maxT} object
#' 
#' @param i \link[base]{logical} or \link[base]{integer} \link[base]{vector},
#' indices of a subset of hypothesis \eqn{\{i_1,\cdots,i_n\}\subset\{1,\cdots,m\}}
#' 
#' @details
#' Function `[.maxT` performs Westfall & Young's \linkS4class{maxT} algorithm  
#' on a subset of test statistics \eqn{\{t_{i_1},\cdots t_{i_n}\}\subset\{t_1,\cdots t_m\}}
#' and their corresponding test statistics 
#' \eqn{\{t_{i_1,b},\cdots,t_{i_n,b}\}\subset\{t_{1,b},\cdots,t_{m,b}\}}
#' in each permuted copy \eqn{b}, \eqn{b=1,\cdots,B}.
#' 
#' In other words, we take the subset of hypothesis \eqn{j=1,\cdots,m},
#' **not** of permutations copies \eqn{b=1,\cdots,B}.
#' 
#' @returns 
#' Function `[.maxT` returns a \linkS4class{maxT} object.
#' 
#' @importFrom methods new
#' @export
'[.maxT' <- function(x, i) {
  new(Class = 'maxT', t. = x@t.[i], T. = x@T.[i, , drop = FALSE], two.sided = x@two.sided)
}


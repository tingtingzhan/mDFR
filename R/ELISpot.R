



#' @title Create `elispot` Object
#' 
#' @description ..
#' 
#' @param formula \link[stats]{formula}
#' 
#' @param data \link[base]{data.frame} with columns in the order of 
#' patient_id, day, antigen, rep1, rep2, etc.
#' 
#' @param control \link[base]{character} scalar, name of negative control wells.
#' Currently do not allow `missing` (for no negative control wells)
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @returns 
#' 
#' Functions [moodie2ELISpot] and [santos2ELISpot] both return an `elispot`.
#' 
#' @examples 
#' # see ?maxT_moodie
#' @references 
#' \url{https://rundfr.fredhutch.org}
#' 
#' @name ELISpot
#' @export
moodie2ELISpot <- function(formula, data, control, ...) {
  if (!is.data.frame(data)) stop('`data` must be data.frame')
  data <- as.data.frame(data) # use S3, 'tibble'
  
  # Remove all-empty response rows # do not do this.  Change DFR function later
  # data <- data[!rowAlls(is.na.data.frame(data[-id.vars]), na.rm = FALSE), , drop = FALSE] 
  # @importFrom matrixStats rowAlls
  
  f_design <- formula[[3L]] # antigen | day + id
  
  xvar <- all.vars(f_design)
  if (anyNA(data[xvar], recursive = TRUE)) stop('Do not allow missingness in ', sQuote(deparse1(f_design)))
  
  f_sort <- call('~', call('+', f_design[[3L]], f_design[[2L]])) # ~Hr + Subj + antigen
  # only symbol `f_design[[2L]]` supported, for now!!
  data <- sort_by.data.frame(data, y = eval(f_sort))
  
  yid <- grepl(pattern = deparse1(formula[[2L]]), x = names(data))
  y <- as.matrix.data.frame(data[yid]) # responses, 'treatment' and 'control'
  dimnames(y) <- NULL
  if (any(y < 0, na.rm = TRUE)) stop('negative value in experiment/control?')
  
  if (!is.character(control) || length(control) != 1L || is.na(control) || !nzchar(control)) stop('illegal `control`')
  
  # `cid`: row indices of negative control
  if (!any(cid <- (data[[f_design[[2L]]]] == control))) stop('Negative control ', sQuote(control), ' is not present in 3rd column of input data')
  
  ret <- data[all.vars(f_design)][!cid, , drop = FALSE] # design *after removing control*
  #X <- data[!cid, id.vars, drop = FALSE] # design matrix *after removing control*
  .rowNamesDF(ret) <- NULL
  
  y1 <- y[!cid, , drop = FALSE]
  dy1 <- dim(y1)
  if (any(id <- (.colSums(is.na(y1), m = dy1[1L], n = dy1[2L]) == dy1[1L]))) {
    y1 <- y1[, !id, drop = FALSE] # remove all-missing columns in `y`
  }
  ret$y1 <- y1
  
  ret$y0 <- y[rep(which(cid), times = tabulate(do.call(interaction, ret[all.vars(f_design[[3L]])]))), ]
  
  attr(ret, which = 'design') <- f_design
  
  class(ret) <- c('elispot', class(ret))
  return(ret)
}






# @param col_trt 'character' scalar (names(data)[3L]), the column name of treatment.  Currently use scalar 'Treatment'
# 
#' @rdname ELISpot
#' 
#' @param ptn0 \link[base]{regex}, pattern of control
#' 
#' @param ptn1 \link[base]{regex}, pattern of treatment
#' 
#' @examples
#' # these are `santos1` and `santos2` in this package
#' santos2ELISpot(formula = ~ antigen | Hr + Subj, data = santos1_raw)
#' santos2ELISpot(formula = ~ antigen | Hr + Subj, data = santos2_raw)
#' 
#' @importFrom zoo na.locf
#' @export
santos2ELISpot <- function(
    formula,
    data, 
    ptn0 = '^y0_',
    ptn1 = '^y1_',
    ...
) {
  if (!is.data.frame(data)) stop('`data` must be a data.frame')
  data <- as.data.frame(data) # use S3, 'tibble'
  
  f_design <- formula[[2L]] # antigen | Hr + Subj
  data[all.vars(f_design)] <- na.locf(data[all.vars(f_design)]) # ?zoo:::na.locf.data.frame
  
  f_sort <- call('~', call('+', f_design[[3L]], f_design[[2L]])) # ~Hr + Subj + antigen
  # only symbol `f_design[[2L]]` supported, for now!!
  
  ret <- sort_by.data.frame(data, y = eval(f_sort))
  .rowNamesDF(ret) <- NULL
  
  nm <- names(ret)
  nm0 <- grepl(pattern = ptn0, x = nm)
  nm1 <- grepl(pattern = ptn1, x = nm)
  y0 <- as.matrix.data.frame(unname(ret[nm0]))
  y1 <- as.matrix.data.frame(unname(ret[nm1]))
  ret[nm0 | nm1] <- NULL
  ret$y1 <- y1
  ret$y0 <- y0
  attr(ret, which = 'design') <- f_design
  class(ret) <- c('elispot', class(ret))
  return(ret)
}





#' @title Several S3 method dispatches for `elispot`
#' 
#' @description
#' `elispot` dispatches for S3 generics \link[base]{log}, \link[base]{split}.
#' 
#' @param x an `elispot`
#' 
#' @param base \link[base]{numeric} scalar, see \link[base]{log}
#' 
# @param f see \link[base]{split.data.frame}, default 
#' 
# @param ... additional parameters, currently not in use
#' 
#' @note
#' Unfortunately \link[base]{log10} and \link[base]{log1p} are not an S3 generics.
#' 
#' \link[base]{log1p} does not have `base` parameter.
#' 
#' @returns
#' 
#' Function [log.elispot] returns an `elispot`.
#' 
#' @name elispot_S3
#' @export log.elispot
#' @export
log.elispot <- function(x, base = exp(1)) {
  if (any(x$y1 == 0, x$y0 == 0, na.rm = TRUE)) {
    x$y1 <- x$y1 + 1
    x$y0 <- x$y0 + 1
  }
  x$y0 <- log(x$y0, base = base)
  x$y1 <- log(x$y1, base = base)
  return(x)
}


# @rdname elispot_S3
# @export split.elispot
# @export
#split.elispot <- function(x, f = eval(call('~', attr(x, which = 'design')[[3L]])), ...) {
#  split.data.frame(x, f = f) # , ...
#  # do not include `...` !! saves a lot of head ache :)
#}




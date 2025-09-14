
#' @title S4 class \linkS4class{elispot}
#' 
#' @slot design \link[base]{data.frame}, experimental design
#' 
#' @slot x1 \link[base]{matrix}
#' 
#' @slot x0 \link[base]{matrix}
#' 
#' @keywords internal
#' @name elispot
#' @aliases elispot-class
#' @export
setClass(Class = 'elispot', slots = c(
  design = 'data.frame',
  x1 = 'matrix',
  x0 = 'matrix'
))



setMethod(f = initialize, signature = 'elispot', definition = function(.Object, ...) {
  
  x <- callNextMethod(.Object, ...)
  
  # remove all-NA columns from `@x0` and `@x1`
  id1 <- (colSums(!is.na(x@x1)) > 0L)
  if (!all(id1)) {
    'all-NA columns in @x1 removed' |>
      message()
    x@x1 <- x@x1[, id1, drop = FALSE]
  }
  
  id0 <- (colSums(!is.na(x@x0)) > 0L)
  if (!all(id0)) {
    'all-NA columns in @x0 removed' |>
      message()
    x@x0 <- x@x0[, , drop = FALSE]
  }

  return(x)
  
})







#' @title Create `elispot` Object, Moodie's fashion
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
#' Functions [moodie2ELISpot()] and [santos2ELISpot()] both return an `elispot`.
#' 
#' @examples 
#' # see ?maxT_moodie
#' @references 
#' \url{https://rundfr.fredhutch.org}
#' 
#' @keywords internal
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
  
  f_sort <- call(name = '~', call(name = '+', f_design[[3L]], f_design[[2L]])) # ~Hr + Subj + antigen
  # only symbol-`f_design[[2L]]` supported, for now!!
  data <- sort_by.data.frame(data, y = eval(f_sort))
  
  xid <- grepl(pattern = deparse1(formula[[2L]]), x = names(data))
  x <- as.matrix.data.frame(data[xid]) # responses, 'treatment' and 'control'
  dimnames(x) <- NULL
  if (any(x < 0, na.rm = TRUE)) stop('negative value in experiment/control?')
  
  if (!is.character(control) || length(control) != 1L || is.na(control) || !nzchar(control)) stop('illegal `control`')
  
  # `cid`: row indices of negative control
  if (!any(cid <- (data[[f_design[[2L]]]] == control))) stop('Negative control ', sQuote(control), ' is not present in 3rd column of input data')
  
  ret <- data[all.vars(f_design)][!cid, , drop = FALSE] # design *after removing control*
  #X <- data[!cid, id.vars, drop = FALSE] # design matrix *after removing control*
  .rowNamesDF(ret) <- NULL
  
  x1 <- x[!cid, , drop = FALSE]
  dx1 <- dim(x1)
  if (any(id <- (.colSums(is.na(x1), m = dx1[1L], n = dx1[2L]) == dx1[1L]))) {
    x1 <- x1[, !id, drop = FALSE] # remove all-missing columns in `x`
  }
  
  x0 <- x[rep(which(cid), times = tabulate(do.call(interaction, ret[all.vars(f_design[[3L]])]))), ]
  
  new(Class = 'elispot', design = ret, x1 = x1, x0 = x0)

}





#' @title Create `elispot` Object, Santos' fashion
#' 
#' @description ..
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
#' @keywords internal
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
  data[all.vars(f_design)] <- data[all.vars(f_design)] |> na.locf() # ?zoo:::na.locf.data.frame
  
  f_sort <- call('~', call('+', f_design[[3L]], f_design[[2L]])) # ~Hr + Subj + antigen
  # only symbol `f_design[[2L]]` supported, for now!!
  
  ret <- sort_by.data.frame(data, y = eval(f_sort))
  .rowNamesDF(ret) <- NULL
  
  nm <- names(ret)
  nm0 <- grepl(pattern = ptn0, x = nm)
  nm1 <- grepl(pattern = ptn1, x = nm)
  x0 <- ret[nm0] |> unname() |> as.matrix.data.frame()
  x1 <- ret[nm1] |> unname() |> as.matrix.data.frame()
  ret[nm0 | nm1] <- NULL
  
  new(Class = 'elispot', design = ret, x1 = x1, x0 = x0)
  
}







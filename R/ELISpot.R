
#' @title S4 class \linkS4class{ELISpot}
#' 
#' @slot design \link[base]{data.frame}, experimental design
#' 
#' @slot x1 \link[base]{matrix}
#' 
#' @slot x0 \link[base]{matrix}
#' 
#' @keywords internal
#' @name ELISpot
#' @aliases ELISpot-class
#' @export
setClass(Class = 'ELISpot', slots = c(
  design = 'data.frame',
  x1 = 'matrix',
  x0 = 'matrix'
))



setMethod(f = initialize, signature = 'ELISpot', definition = function(.Object, ...) {
  
  x <- callNextMethod(.Object, ...)
  
  # remove all-NA columns from `@x0` and `@x1`
  .rm_colNA <- \(x) { # `x` is 'matrix'
    id <- (colSums(!is.na(x)) > 0L)
    if (all(id)) return(x)
    substitute(x) |>
      deparse1() |>
      col_red() |>
      sprintf(fmt = 'all-NA columns in %s removed') |>
      message()
    x <- x[, id, drop = FALSE]
    return(x)
  }
   
  x@x1 <- .rm_colNA(x@x1)
  x@x0 <- .rm_colNA(x@x0)

  return(x)
  
})







#' @title Create `ELISpot` Object, Moodie's fashion
#' 
#' @description ..
#' 
#' @param formula \link[stats]{formula}
#' 
#' @param data \link[base]{data.frame} with columns in the order of 
#' patient_id, day, antigen, rep1, rep2, etc.
#' 
#' @param pattern1 \link[base]{regex}
#' 
#' @param pattern0 \link[base]{character} scalar, name of negative control wells.
#' Currently do not allow `missing` (for no negative control wells)
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @returns 
#' 
#' Functions [santos2ELISpot()] both return an `ELISpot`.
#' 
#' @keywords internal
#' @export
as.ELISpot <- function(
    formula, 
    data,
    pattern0,
    pattern1,
    ...
) {
  
  if (!is.data.frame(data)) stop('`data` must be data.frame')
  data <- as.data.frame(data) # use S3, 'tibble'
  
  if (!is.character(pattern0) || length(pattern0) != 1L || is.na(pattern0) || !nzchar(pattern0)) stop('illegal `pattern0`')
  if (!is.character(pattern1) || length(pattern1) != 1L || is.na(pattern1) || !nzchar(pattern1)) stop('illegal `pattern1`')
  
  # NOT remove all-empty response rows here

  f_design <- formula[[2L]] # antigen | day + id
  if (f_design[[1L]] != '|') stop('design formula must contain a `|`')
   
  xvar <- all.vars(f_design)
  if (anyNA(data[xvar], recursive = TRUE)) stop('Do not allow missingness in ', sQuote(deparse1(f_design)))
  
  # do NOT sort by antigen!!
  data <- call(name = '~', call(name = '+', f_design[[3L]], f_design[[2L]])) |> 
    #call(name = '~', f_design[[3L]]) |> 
    # ~Hr + Subj + antigen
    # only symbol-`f_design[[2L]]` supported, for now!!
    eval() |>
    sort_by.data.frame(x = data, y = _)

  X <- data |>
    names() |>
    grepv(pattern = pattern1, x = _) |>
    subset.data.frame(x = data, select = _) |>
    as.matrix.data.frame()
  # responses, 'treatment'-only, or 'treatment'-'control'
  # was named `x`
  dimnames(X) <- NULL
  if (any(X < 0, na.rm = TRUE)) stop('negative value in experiment/control?')
  
  v0 <- data |>
    names() |>
    grepv(pattern = pattern0, x = _)
  
  if (!length(v0)) { # moodie-style
    
    # `cid`: row indices of negative control
    if (!any(cid <- (data[[f_design[[2L]]]] == pattern0))) stop('Negative control ', sQuote(pattern0), ' is not present in data_design')
    
    ret <- data[all.vars(f_design)][!cid, , drop = FALSE] # design *after removing control*
    .rowNamesDF(ret) <- NULL
    
    x1 <- X[!cid, , drop = FALSE]
    dx1 <- dim(x1)
    if (any(id <- (.colSums(is.na(x1), m = dx1[1L], n = dx1[2L]) == dx1[1L]))) {
      x1 <- x1[, !id, drop = FALSE] # remove all-missing columns in `x`
    }
    
    x0 <- X[rep(which(cid), times = tabulate(do.call(interaction, ret[all.vars(f_design[[3L]])]))), ]
    
    return(new(Class = 'ELISpot', design = ret, x1 = x1, x0 = x0))
    
  }
  

}





#' @title Create `ELISpot` Object, Santos' fashion
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
  
  new(Class = 'ELISpot', design = ret, x1 = x1, x0 = x0)
  
}







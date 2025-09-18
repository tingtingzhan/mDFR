
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

  # Both `@x1` and `@x0` must contain more than 1 measurement per row
  id <- (rowSums(!is.na(x@x1)) > 1L) & (rowSums(!is.na(x@x0)) > 1L)
  if (!all(id)) {
    x@design <- x@design[id, ]
    x@x1 <- x@x1[id, , drop = FALSE]
    x@x0 <- x@x0[id, , drop = FALSE]
  }
  
  return(x)
  
})






#' @title Create \linkS4class{ELISpot} Object
#' 
#' @description ..
#' 
#' @param formula \link[stats]{formula}
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param pattern1 \link[base]{regex}
#' 
#' @param pattern0 \link[base]{character} scalar of name of negative-control arm, or of \link[base]{regex} control pattern
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @returns 
#' Functions [as.ELISpot()] both return an \linkS4class{ELISpot} object
#' 
#' @keywords internal
#' @importFrom zoo na.locf
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
  data[xvar] <- data[xvar] |> 
    na.locf() # ?zoo:::na.locf.data.frame
  if (anyNA(data[xvar], recursive = TRUE)) stop('Do not allow missingness in ', sQuote(deparse1(f_design)))
  
  # NEXT: do NOT sort by antigen!!
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
  
  if (!length(v0)) { # Moodie-style
    
    # `cid`: row indices of negative control
    if (!any(cid <- (data[[f_design[[2L]]]] == pattern0))) stop('Negative control ', sQuote(pattern0), ' is not present in data_design')
    
    ret <- data[xvar][!cid, , drop = FALSE] # design *after removing control*
    .rowNamesDF(ret) <- NULL
    
    x1 <- X[!cid, , drop = FALSE]
    dx1 <- dim(x1)
    if (any(id <- (.colSums(is.na(x1), m = dx1[1L], n = dx1[2L]) == dx1[1L]))) {
      x1 <- x1[, !id, drop = FALSE] # remove all-missing columns in `x`
    }
    
    x0 <- X[rep(which(cid), times = tabulate(do.call(interaction, ret[all.vars(f_design[[3L]])]))), ]
    
    return(new(Class = 'ELISpot', design = ret, x1 = x1, x0 = x0))
    
  } # Moodie-style
  
  # Santos-style
  
  x0 <- v0 |>
    subset.data.frame(x = data, select = _) |>
    as.matrix.data.frame()
  
  return(new(Class = 'ELISpot', design = data[xvar], x1 = X, x0 = x0))
  
}



#' @title Permutation Indices of Treatment in \linkS4class{ELISpot}
#' 
#' @description
#' Permuted indices of treatment in an \linkS4class{ELISpot} object
#' 
#' @note
#' Westfall & Young's \linkS4class{maxT} algorithm
#' (Box 2, page 82 of \doi{10.1214/ss/1056397487}) says *permutation*.
#' 
#' In **R** convention, this is a *combination* \link[utils]{combn}.
#' 
#' @param x an \linkS4class{ELISpot} object
#' 
#' @returns 
#' Function [combn_ELISpot()] returns a \link[base]{list} of \link[base]{integer} \link[base]{vector}s.
#' 
#' @keywords internal
#' @importFrom utils combn
#' @export
combn_ELISpot <- function(x) {
  n1 <- dim(x@x1)[2L]
  n0 <- dim(x@x0)[2L]
  combn(n1 + n0, m = n1, simplify = FALSE)
}






#' @title \link[base]{log} of \linkS4class{ELISpot}
#' 
#' @description
#' \linkS4class{ELISpot} dispatches for S3 generics \link[base]{log}.
#' 
#' @param x an \linkS4class{ELISpot} object
#' 
#' @param base \link[base]{numeric} scalar, see \link[base]{log}
#' 
#' @note
#' The `S3` generic function \link[base]{log1p} does not have `base` parameter.
#' 
#' @returns
#' 
#' Function [log.ELISpot()] returns an `ELISpot`.
#' 
#' @keywords internal
#' @name log_ELISpot
#' @export log.ELISpot
#' @export
log.ELISpot <- function(x, base = exp(1)) {
  if (any(x@x1 == 0, x@x0 == 0, na.rm = TRUE)) {
    x@x1 <- x@x1 + 1
    x@x0 <- x@x0 + 1
  }
  x@x0 <- log(x@x0, base = base)
  x@x1 <- log(x@x1, base = base)
  return(x)
}


#' @rdname log_ELISpot
#' @export log1p.ELISpot
#' @export
log1p.ELISpot <- function(x) {
  x@x0 <- log1p(x@x0)
  x@x1 <- log1p(x@x1)
  return(x)
}

#' @rdname log_ELISpot
#' @export log10.ELISpot
#' @export
log10.ELISpot <- function(x) {
  x@x0 <- log10(x@x0)
  x@x1 <- log10(x@x1)
  return(x)
}

#' @rdname log_ELISpot
#' @export log2.ELISpot
#' @export
log2.ELISpot <- function(x) {
  x@x0 <- log2(x@x0)
  x@x1 <- log2(x@x1)
  return(x)
}








#' @title \link[base]{split} an \linkS4class{ELISpot} Object
#' 
#' @param x an \linkS4class{ELISpot}
#' 
#' @param f ..
#' 
#' @param drop \link[base]{logical} scalar, must be `FALSE`
#' 
#' @keywords internal
#' @export split.ELISpot
#' @export
split.ELISpot <- function(x, f, drop = TRUE, ...) {
  
  if (!drop) stop('`drop` must be TRUE')
  
  # most simple solution..
  z <- x@design
  z$x1 <- x@x1
  z$x0 <- x@x0
  
  foo <- \(y) {
    y0 <- y
    y0$x1 <- y0$x0 <- NULL
    new(Class = 'ELISpot', design = y0, x1 = y$x1, x0 = y$x0)
  }

  z |>
    split.data.frame(f = f, drop = drop, ...) |>
    lapply(FUN = foo)
  
}






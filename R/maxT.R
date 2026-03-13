
#' @title \linkS4class{maxT} of Distribution Free Statistics
#' 
#' @description ..
#' 
#' @param x see **Usage**
#' 
# @param alternative \link[base]{logical} scalar, see \linkS4class{maxT}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @references 
#' \url{https://tingtingzhan-maxt.netlify.app/appendix/westfall_yang.html}
#' 
#' @keywords internal
#' @name maxT
#' @export
maxT <- function(x, ...) UseMethod(generic = 'maxT')


#' @rdname maxT
#' @export
maxT.free_d <- function(x, ...) {
  
  x@data |> # supported: 'ELISpot'
    permID() |>
    lapply(FUN = free_d, x = x@data, s4 = FALSE) |>
    do.call(what = cbind, args = _) |>
    new(
      Class = 'maxT', 
      t. = x@.Data, T. = _,
      design = x@data@design,
      label = labels(x@data), # [labels.ELISpot()], etc
      ...
    )

}




#' @rdname maxT
#' @export
maxT.free_t <- function(x, ...) {
  
  x@data |>
    permID() |>
    lapply(FUN = free_t, x = x@data, s4 = FALSE) |>
    do.call(what = cbind, args = _) |>
    new(
      Class = 'maxT', 
      t. = x@.Data, T. = _,
      design = x@data@design,
      label = labels(x@data), # [labels.ELISpot()], etc
      ...
    )

}





#' @rdname maxT
#' @importFrom matrixStats colAnys
#' @export
maxT.free_t_diff <- function(x, ...) {
  
  nr <- nrow(x@e1@data@design) # `nrow(x@e1@data) == nrow(data0)` enforced in [`-`('free_t', 'free_t')]
  
  t_ <- free_t(x@e1@data) - free_t(x@e2@data)

  # based on permutation
  e1 <- x@e1@data |>
    permID() |>
    lapply(FUN = free_t, x = x@e1@data, s4 = FALSE)
  n1 <- length(e1)
  e2 <- x@e2@data |> 
    permID() |>
    lapply(FUN = free_t, x = x@e2@data, s4 = FALSE)
  n0 <- length(e2)
  
  ### actually [`-`('free_t', 'free_t')]
  d1. <- e1 |>
    lapply(FUN = attr, which = 'delta', exact = TRUE) |>
    unlist(use.names = FALSE)
  df1. <- e1 |>
    lapply(FUN = attr, which = 'df', exact = TRUE) |>
    unlist(use.names = FALSE)
  stderr1. <- e1 |>
    lapply(FUN = attr, which = 'stderr', exact = TRUE) |>
    unlist(use.names = FALSE)
  dim(d1.) <- dim(df1.) <- dim(stderr1.) <- c(nr, n1)
  
  d2. <- e2 |> 
    lapply(FUN = attr, which = 'delta', exact = TRUE) |>
    unlist(use.names = FALSE)
  df2. <- e2 |>
    lapply(FUN = attr, which = 'df', exact = TRUE) |>
    unlist(use.names = FALSE)
  stderr2. <- e2 |>
    lapply(FUN = attr, which = 'stderr', exact = TRUE) |>
    unlist(use.names = FALSE)
  dim(d2.) <- dim(df2.) <- dim(stderr2.) <- c(nr, n0)
  
  id_ <- expand.grid(e2 = seq_len(n0), e1 = seq_len(n1))
  var_pooled <- ((df1.*stderr1.^2)[,id_$e1] + (df2.*stderr2.^2)[,id_$e2]) / (df1.[,id_$e1] + df2.[,id_$e2])
  T. <- (d1.[,id_$e1] - d2.[,id_$e2]) / sqrt(var_pooled)
  if (any(id <- colAnys(is.na(T.)))) {
    message(sprintf(fmt = '%.1f%% permutations with NA test-statistics are omitted', 1e2*mean.default(id)))
    T. <- T.[, !id]
  }
  ### but [`-`('free_t', 'free_t')] is too slow on \emph{permutation-of-permutation}
  ### have to manually vectorize !!
  
  new(
    Class = 'maxT', 
    t. = t_, T. = T.,
    design = x@design,
    label = labels(x), # [labels.free_t_diff()]
    ...
  )
  
}




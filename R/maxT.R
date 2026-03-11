
#' @title \linkS4class{maxT} of Distribution Free Statistics
#' 
#' @description ..
#' 
#' @param x see **Usage**
#' 
#' @param two.sided \link[base]{logical} scalar, see \linkS4class{maxT}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @references 
#' \url{https://tingtingzhan-maxt.netlify.app/appendix/westfall_yang.html}
#' 
#' @keywords internal
#' @name maxT
#' @export
maxT <- function(x, two.sided, ...) UseMethod(generic = 'maxT')


#' @rdname maxT
#' @export
maxT.free_d <- function(x, two.sided = TRUE, ...) {
  
  x@data |> # supported: 'ELISpot'
    permID() |>
    lapply(FUN = free_d, x = x@data, s4 = FALSE) |>
    do.call(what = cbind, args = _) |>
    new(
      Class = 'maxT', 
      t. = x@.Data, T. = _,
      design = x@data@design,
      name = labels(x@data),
      two.sided = two.sided
    )

}




#' @rdname maxT
#' @export
maxT.free_t <- function(x, two.sided = TRUE, ...) {
  
  x@data |>
    permID() |>
    lapply(FUN = free_t, x = x@data, s4 = FALSE) |>
    do.call(what = cbind, args = _) |>
    new(
      Class = 'maxT', 
      t. = x@.Data, T. = _,
      design = x@data@design,
      name = labels(x@data),
      two.sided = two.sided
    )

}





#' @rdname maxT
#' @importFrom matrixStats colAnys
#' @export
maxT.free_t_diff <- function(x, two.sided = TRUE, ...) {
  
  nr <- nrow(x@e1@data@design) # `nrow(x@e1@data) == nrow(data0)` enforced in [`-`('free_t', 'free_t')]
  
  t_ <- free_t(x@e1@data) - free_t(x@e2@data)

  # based on permutation
  tm1 <- x@e1@data |>
    permID() |>
    lapply(FUN = free_t, x = x@e1@data, s4 = FALSE)
  n1 <- length(tm1)
  tm0 <- x@e2@data |> 
    permID() |>
    lapply(FUN = free_t, x = x@e2@data, s4 = FALSE)
  n0 <- length(tm0)
  
  ### actually [`-`('free_t', 'free_t')]
  d1. <- tm1 |>
    lapply(FUN = attr, which = 'delta', exact = TRUE) |>
    unlist(use.names = FALSE)
  df1. <- tm1 |>
    lapply(FUN = attr, which = 'df', exact = TRUE) |>
    unlist(use.names = FALSE)
  stderr1. <- tm1 |>
    lapply(FUN = attr, which = 'stderr', exact = TRUE) |>
    unlist(use.names = FALSE)
  dim(d1.) <- dim(df1.) <- dim(stderr1.) <- c(nr, n1)
  
  d0. <- tm0 |> 
    lapply(FUN = attr, which = 'delta', exact = TRUE) |>
    unlist(use.names = FALSE)
  df0. <- tm0 |>
    lapply(FUN = attr, which = 'df', exact = TRUE) |>
    unlist(use.names = FALSE)
  stderr0. <- tm0 |>
    lapply(FUN = attr, which = 'stderr', exact = TRUE) |>
    unlist(use.names = FALSE)
  dim(d0.) <- dim(df0.) <- dim(stderr0.) <- c(nr, n0)
  
  id_ <- expand.grid(tm0 = seq_len(n0), tm1 = seq_len(n1))
  var_pooled <- ((df1.*stderr1.^2)[,id_$tm1] + (df0.*stderr0.^2)[,id_$tm0]) / (df1.[,id_$tm1] + df0.[,id_$tm0])
  T. <- (d1.[,id_$tm1] - d0.[,id_$tm0]) / sqrt(var_pooled)
  if (any(id <- colAnys(is.na(T.)))) {
    message(sprintf(fmt = '%.1f%% permutations with NA test-statistics are omitted', 1e2*mean.default(id)))
    T. <- T.[, !id]
  }
  ### but [`-`('free_t', 'free_t')] is too slow on \emph{permutation-of-permutation}
  ### have to manually vectorize !!
  
  new_design <- mapply(
    FUN = \(c1, c0) { # operation per-*c*olumn
      if (anyNA(c1) || anyNA(c0)) stop('does not allow NA in experiment design')
      if (all(c1 == c0)) return(c1)
      return(paste(c1, c0, sep = ' vs. '))
    }, 
    c1 = x@e1@data@design, 
    c0 = x@e2@data@design, 
    SIMPLIFY = FALSE
  ) |>
    as.data.frame.list()
  
  new(
    Class = 'maxT', 
    t. = t_, T. = T.,
    design = new_design,
    name = labels_design(new_design),
    two.sided = two.sided
  )
  
}





#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param x an \linkS4class{ELISpot}
#' 
#' @param two.sided \link[base]{logical} scalar, see \linkS4class{maxT}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @references 
#' 
#' Radleigh Santos's 2015 paper \doi{10.3390/cells4010001}
#' 
#' @keywords internal
#' @name maxT
#' @export
maxT <- function(x, two.sided, ...) UseMethod(generic = 'maxT')


#' @rdname maxT
#' @export maxT.free_d
#' @export
maxT.free_d <- function(x, two.sided = TRUE, ...) {
  
  T_ <- x@data |>
    combn_ELISpot() |>
    lapply(FUN = free_d.ELISpot, x = x@data, s4 = FALSE) |>
    unlist(use.names = FALSE)
  nr <- nrow(x@data@design)
  dim(T_) <- c(nr, length(T_)/nr)
  
  tmp <- x@data@design |>
    lapply(FUN = \(i) {
      if (is.factor(i)) i <- as.character.factor(i)
      unique(i)
    })
  
  new(
    Class = 'maxT', 
    t. = x@.Data, T. = T_,
    design = x@data@design,
    name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; '),
    two.sided = two.sided
  )
}





#' @rdname maxT
#' @export maxT.free_t
#' @export
maxT.free_t <- function(x, two.sided = TRUE, ...) {
  
  T_ <- x@data |>
    combn_ELISpot() |>
    lapply(FUN = free_t.ELISpot, x = x@data, s4 = FALSE) |>
    unlist(use.names = FALSE)
  nr <- nrow(x@data@design)
  dim(T_) <- c(nr, length(T_)/nr)
  
  tmp <- x@data@design |>
    lapply(FUN = \(i) {
      if (is.factor(i)) i <- as.character.factor(i)
      unique(i)
    })
  
  new(
    Class = 'maxT', 
    t. = x@.Data, T. = T_,
    design = x@data@design,
    name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; '),
    two.sided = two.sided
  )
}







# @title Modified Distribution Free Resampling, at Two Timepoints
# @param ... additional parameters, such as `null.value` in function [free_t] and `two.sided` for \linkS4class{maxT}

#' @rdname maxT
#' @importFrom matrixStats colAnys
#' @export maxT.free_t_diff
#' @export
maxT.free_t_diff <- function(x, two.sided = TRUE, ...) {
  
  if (nrow(x@e1@data@design) != nrow(x@e2@data@design)) {
    # timepoint1 and timepoint2 may not have same subjects!!
    
    # before 2025-09
    #d_1 <- x@e1@data; d_1$x1 <- d_1$x0 <- NULL
    #d_0 <- x@e2@data; d_0$x1 <- d_0$x0 <- NULL
    #tmp <- mapply(FUN = intersect, d_1, d_0)
    #d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    #data1 <- merge.data.frame(d_share, x@e1@data, by = names(d_share), all.x = TRUE)
    #data0 <- merge.data.frame(d_share, x@e2@data, by = names(d_share), all.x = TRUE)
    # end of before 2025-09
    
    #d_1 <- x@e1@data@design
    #d_0 <- x@e2@data@design
    #tmp <- mapply(FUN = intersect, d_1, d_0)
    #d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    stop('um, tzh needs to think about how to write this beautifully..')
    #data1 <- merge.data.frame(d_share, x@e1@data, by = names(d_share), all.x = TRUE)
    #data0 <- merge.data.frame(d_share, x@e2@data, by = names(d_share), all.x = TRUE)
    
  }
  
  nr <- nrow(x@e1@data@design) # now `nrow(x@e1@data) == nrow(data0)`
  
  t_ <- free_t.ELISpot(x@e1@data) - free_t.ELISpot(x@e2@data)

  # based on permutation
  tm1 <- x@e1@data |>
    combn_ELISpot() |>
    lapply(FUN = free_t.ELISpot, x = x@e1@data, s4 = FALSE)
  n1 <- length(tm1)
  tm0 <- x@e2@data |> 
    combn_ELISpot() |>
    lapply(FUN = free_t.ELISpot, x = x@e2@data, s4 = FALSE)
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
  
  # combine `x@e1@data` and `x@e2@data` for output
  d1 <- x@e1@data@design
  d0 <- x@e2@data@design
  if (!identical(names(d1), names(d0))) stop('`ELISpot` at two time points must have same design')
  d <- mapply(FUN = \(c1, c0) {
    if (anyNA(c1) || anyNA(c0)) stop('does not allow NA in experiment design')
    if (all(c1 == c0)) return(c1)
    return(paste(c1, c0, sep = ' vs. '))
  }, c1 = d1, c0 = d0, SIMPLIFY = FALSE)
  tmp <- lapply(d, FUN = unique.default)
  
  new(
    Class = 'maxT', 
    t. = t_, T. = T.,
    design = as.data.frame.list(d),
    name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; '),
    two.sided = two.sided
  )
  
}

























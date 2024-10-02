

trim_elispot <- function(data) {
  id <- (rowSums(!is.na(data$y1)) > 1L) & (rowSums(!is.na(data$y1)) > 1L)
  return(data[id, ])
}
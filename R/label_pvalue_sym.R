

#' Function [label_pvalue_sym()] is used so extensively, 
#' therefore tzh keeps it in \pkg{rmd.tzh} instead of \pkg{scales.tzh}.



#' @title Label \eqn{p}-values with Significance Symbol
#' 
#' @description
#' Label \eqn{p}-values with significance symbol.
#' 
#' @param ... parameters of function \link[scales]{label_pvalue}
#' 
#' @details
#' 
#' Pipeline:
#' 
#' Step 1: Apply function \link[scales]{label_pvalue}.
#' 
#' Step 2: Drop leading zeros.
#' 
#' Step 3: Apply function \link[stats]{symnum} (also see function \link[stats]{printCoefmat}).
#' 
#' @note
#' Function \link[scales]{label_pvalue} is much prettier and more flexible than function \link[base]{format.pval}.
#' 
#' @returns 
#' Function [label_pvalue_sym()] returns a \link[base]{function}.
#' 
#' @examples 
#' p = c(a = pi^-100, b = .02, c = .05, d = .1, e = .9999, f = NA_real_)
#' p |> label_pvalue_sym()()
#' p |> label_pvalue_sym(add_p = TRUE)()
#' 
#' # below: exception handling
#' double() |> scales::label_pvalue()() # do not like!
#' # ?scales::pvalue is bad; ?scales::number is fine!
#' double() |> label_pvalue_sym()() # nice!
#' @keywords internal
#' @importFrom scales label_pvalue
#' @importFrom stats symnum
#' @export
label_pvalue_sym <- function(...) {
  
  function(x) { # see ?scales::label_pvalue; `...` no need to be in the args
    
    ret <- x
    storage.mode(ret) <- 'character'
    # `attributes(x)` kept
    
    if (!length(x)) return(ret)
    
    ret[] <- x |> 
      label_pvalue(...)() |>
      sub(pattern = '([-]?)0[.]', replacement = '\\1.') # http://stackoverflow.com/questions/12643391
    
    if (getOption('show.signif.stars')) { # see ?stats::printCoefmat
      sym <- symnum(
        x, corr = FALSE, na = FALSE, 
        cutpoints = c(0, .001, .01, .05, .1, 1), 
        symbols = 
          # c("***", "**", "*", ".", " ")
          # c('\u2605\u2605\u2605', '\u2605\u2605', '\u2605', '\u2606', '') # star
          c('\u2b51\u2b51\u2b51', '\u2b51\u2b51', '\u2b51', '\u2b52', '') # small star
      ) # see ?stats::printCoefmat
      ret[] <- ret |> 
        paste(sym) |> 
        trimws()
    }
    
    ret[is.na(x)] <- '' # *not* NA_character_
    
    return(ret)
      
  }

}



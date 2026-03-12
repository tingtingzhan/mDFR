

# The function [label_pvalue_sym()] is used so extensively, 
# therefore tzh keeps it in package \pkg{fastmd} instead of package \pkg{scales.tzh}.
 



#' @title Label \eqn{p}-values with Significance Symbol
#' 
#' @description
#' Label \eqn{p}-values with significance symbol.
#' 
#' @param star \link[base]{character} \link[base]{vector} 
#' of \link[base]{length}-`2L` to denote the significance levels. 
#' For example, to use \link[stats]{printCoefmat}-style of significance notation,
#' users should let `star = c('*', '.')`
#' 
#' @param ... parameters of the function \link[scales]{label_pvalue}
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
#' The function \link[scales]{label_pvalue} is much prettier and more flexible than function \link[base]{format.pval}.
#' 
#' @returns 
#' The function [label_pvalue_sym()] returns a \link[base]{function}.
#' 
#' @examples 
#' p = c(a = pi^-100, b = .02, c = .05, d = .1, e = .9999, f = NA_real_)
#' p |> label_pvalue_sym()() # star
#' p |> label_pvalue_sym(star = c('\u2605', '\u2606'))() # small star
#' p |> label_pvalue_sym(star = c('*', '.'))() # ?stats::printCoefmat
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
label_pvalue_sym <- function(star = c('\u2b51', '\u2b52'), ...) {
  
  \(x) { # see ?scales::label_pvalue; parameters no need to be in the args!!
    
    ret <- x
    storage.mode(ret) <- 'character'
    # `attributes(x)` kept
    
    if (!length(x)) return(ret)
    
    ret[] <- x |> 
      label_pvalue(...)() |>
      sub(pattern = '([-]?)0[.]', replacement = '\\1.', x = _) # http://stackoverflow.com/questions/12643391
    
    if (getOption('show.signif.stars')) { # see ?stats::printCoefmat
      sym <- symnum(
        x, corr = FALSE, na = FALSE, 
        cutpoints = c(0, .001, .01, .05, .1, 1), 
        symbols = c( # do not want to import ?stringi::stri_dup
          star[1L] |> rep(times = 3L) |> paste0(collapse = ''),
          star[1L] |> rep(times = 2L) |> paste0(collapse = ''),
          star[1L],
          star[2L],
          ''
        )
      )
      ret[] <- ret |> 
        paste(. = _, sym) |> 
        trimws()
    }
    
    ret[is.na(x)] <- '' # *not* NA_character_
    
    return(ret)
      
  }

}



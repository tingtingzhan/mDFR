% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxT.R
\name{autoplot.maxT}
\alias{autoplot.maxT}
\alias{autoplot_maxT_}
\title{Visualize Westfall & Young's \linkS4class{maxT} Algorithm}
\usage{
\method{autoplot}{maxT}(object, ...)

autoplot_maxT_(
  p_perm,
  p_mono,
  tr,
  U,
  two.sided,
  conf.level = 0.95,
  vertical = TRUE,
  ...
)
}
\arguments{
\item{object}{a \linkS4class{maxT} object}

\item{...}{additional parameters, currently not in use}

\item{p_perm, p_mono, tr, U, two.sided}{description}

\item{conf.level}{\link[base]{double} scalar,
confidence level, or \eqn{1-\alpha}. Default .95 (or \eqn{\alpha=.05})}

\item{vertical}{\link[base]{logical} scalar, whether the
labels of ordered test statistics and permutation adjusted \eqn{p}-values should be on vertical axes,
default \code{TRUE}}
}
\value{
Function \link{autoplot.maxT} returns a \link[ggplot2]{ggplot} object.
}
\description{
To visualize Westfall & Young's \linkS4class{maxT} algorithm using package \CRANpkg{ggplot2}.
}
\details{
Function \link{autoplot.maxT} plots
\itemize{
\item the successive maxima
\eqn{u_{jb}}, \eqn{j=1,\cdots,m}, \eqn{b=1,\cdots,B}, and
\item the decreasing-ordered statistics \eqn{|t_{r_1}|\geq|t_{r_2}|\geq\cdots\geq|t_{r_m}|} for two-sided test,
or \eqn{t_{r_1}\geq t_{r_2}\geq\cdots\geq t_{r_m}} for one-sided test
}
Printed on opposing axis are
\itemize{
\item values of the decreasing-ordered statistics
\item permutation adjusted \eqn{p}-values \eqn{\tilde{p}_{r_j}}, as well as under monotonicity constraints \eqn{\tilde{p}^*_{r_j}}
}
Tests with \eqn{\tilde{p}^*_{r_j}\leq\alpha} is considered significant
and colored pink (hex color \verb{#F8766D}), otherwise non-significant and colored blue (hex color \verb{#00BFC4})

See full details of these notations in \linkS4class{maxT}.
}

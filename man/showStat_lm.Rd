% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/showStat.R
\name{showStat_lm}
\alias{showStat_lm}
\title{Draw distribution of statistic for an object of lm}
\usage{
showStat_lm(fit)
}
\arguments{
\item{fit}{An object of class lm}
}
\description{
Draw distribution of statistic for an object of lm
}
\examples{
fit=lm(mpg~hp*wt+vs,data=mtcars)
showStat_lm(fit)
}

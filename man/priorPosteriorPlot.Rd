% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reclassify.R
\name{priorPosteriorPlot}
\alias{priorPosteriorPlot}
\title{Draw prior versus posterior plot}
\usage{
priorPosteriorPlot(fitlist, xlab, ylab, title, ...)
}
\arguments{
\item{fitlist}{A list of objects of class glm}

\item{xlab, ylab, title}{character}

\item{...}{Further arguments to be passed to geom_point()}
}
\description{
Draw prior versus posterior plot
}
\examples{
form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
priorPosteriorPlot(list(fit1,fit2),alpha=0.3,color="red")
}

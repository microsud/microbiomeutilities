% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkg_utils.R
\name{taxa_pooler_mcola}
\alias{taxa_pooler_mcola}
\title{Pool Taxa}
\usage{
taxa_pooler_mcola(x)
}
\arguments{
\item{x}{\code{\link{phyloseq-class}} object}
}
\description{
Creates a list of dataframes at different taxonomic levels.
}
\examples{
# library(phyloseq)
# library(microbiome)
# data("zackular2014")
# sub.sm <- sample(sample_names(zackular2014), 20)
# pseq <- prune_samples(sub.sm,zackular2014)
# taxa_pooler_mcola(pseq)
}
\keyword{utilities}

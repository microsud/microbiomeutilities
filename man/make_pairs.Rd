% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkg_utils.R
\name{make_pairs}
\alias{make_pairs}
\title{Make pairs}
\usage{
make_pairs(x)
}
\arguments{
\item{x}{list of vector to compare}
}
\description{
Creates a combination of variables
for use with ggpubr::stat_compare_means.
}
\examples{
# library(microbiomeutilities)
# data("zackular2014")
# pseq <- zackular2014
# comps <- make_pairs(meta(pseq)$DiseaseState)
}
\keyword{utilities}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_paired_abundances.R
\name{plot_paired_abundances}
\alias{plot_paired_abundances}
\title{A paired-boxplot for user specified list of taxa}
\usage{
plot_paired_abundances(
  x,
  select.taxa = NULL,
  group = NULL,
  group.colors = NULL,
  dot.opacity = 0.25,
  dot.size = 2,
  add.box = FALSE,
  box.opacity = 0.25,
  group.order = NULL,
  add.violin = TRUE,
  violin.opacity = 0.25,
  ncol = NULL,
  nrow = NULL,
  line = NULL,
  line.down = "#7209b7",
  line.stable = "#8d99ae",
  line.up = "#14213d",
  line.na.value = "grey50",
  line.guide = "legend",
  line.opacity = 0.25,
  line.size = 1,
  jitter.width = 0
)
}
\arguments{
\item{x}{\code{\link{phyloseq-class}} object.}

\item{select.taxa}{a character list of taxa to be plotted. eg. select.taxa <- c("OTU-370251", "OTU-311173", "OTU-341024").}

\item{group}{Grouping variable to compare. x axis, eg. before-after, t1-t2.}

\item{group.colors}{Colors for plotting groups.}

\item{dot.opacity}{For ggplot alpha to determine opacity for points.}

\item{dot.size}{For ggplot point size.}

\item{add.box}{Logical. If boxplot to the added. Default=TRUE}

\item{box.opacity}{For ggplot alpha to determine opacity for box.}

\item{group.order}{Default is NULL. a list specifying order of x-axis.}

\item{add.violin}{Logical. If half violin to the added. Default=TRUE}

\item{violin.opacity}{If add.violin=TRUE, opacity for violin.}

\item{ncol}{If 2 or more taxa to plot, specify number of columns.}

\item{nrow}{If 2 or more taxa to plot, specify number of rows.}

\item{line}{Variable to use for lines. E.g. "subject" before-after}

\item{line.down}{Line Color for change when negative. Decreased abundance.}

\item{line.stable}{Line Color for no change.}

\item{line.up}{Line Color for change when positive. Increased abundance.}

\item{line.na.value}{"grey50" for no/missing observations.}

\item{line.guide}{"none" to plot guide for line.}

\item{line.opacity}{Line opacity.}

\item{line.size}{Size of line to plot.}

\item{jitter.width}{Value to avoid over plotting by moving points.}
}
\value{
\code{\link{ggplot}} object. This can be further modified using ggpubr.
}
\description{
User specified taxa are plotted.
}
\details{
Useful for instances where user is interested only in some taxa and thier change after
         an intervention. This can also be used at higher taxonomic
         levels, output from phyloseq::tax_glom or microbiome::aggregate_taxa.
}
\examples{
\dontrun{
library(microbiome)
library(microbiomeutilities)
library(gghalves)
library(tidyr)
data(peerj32) # Source: https://peerj.com/articles/32/
pseq <- peerj32$phyloseq # Ren
pseq.rel <- microbiome::transform(pseq, "compositional")
select.taxa <- c("Akkermansia", "Dialister")
group.colors <- c("brown3", "steelblue", "grey70")
p <- plot_paired_abundances(pseq.rel,
  select.taxa = select.taxa,
  group = "time",
  group.colors = group.colors,
  dot.opacity = 0.25,
  dot.size = 2,
  group.order = NULL,
  line = "subject"
)
p
}
}
\keyword{visualization}

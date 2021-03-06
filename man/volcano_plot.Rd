% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick_plots.R
\name{volcano_plot}
\alias{volcano_plot}
\title{Quickly plot a volcano plot}
\usage{
volcano_plot(fit, fc_cutoff = 2, pval_cutoff = 0.05)
}
\arguments{
\item{fit}{tidied dataframe of fit results}

\item{fc_cutoff}{cutoff to determine what is considered differentially expressed}

\item{pval_cutoff}{plots horizontal line at p-value}
}
\value{
a `list` containing a volcano plot (ggplot object) and the data used to
generate the volcano plot
}
\description{
`volcano_plot()` plots a volcano plot
}
\details{
Like any other ggplot object, you can customize the theme of the plot.
Note that this function generates a fairly generic volcano plot, and thus is
intended for quick exploratory purposes (much like the intention behind `qplot()`).
}
\examples{
# As an EDA step within a pipeline of functions

counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") \%>\%
  mutate(FileName = stringr::str_replace(FileName, "\\\\.", "-"))

id <- as.character(counts$EntrezGeneID)
my_design <- check_sample_names(counts, c(1,2), meta, FileName) \%>\%
  purrr::pluck("meta") \%>\%
  make_design_matrix(., c("Status"))

check_sample_names(counts, c(1,2), meta, FileName) \%>\%
  purrr::pluck("mod_count") \%>\%
  filter_genes(., id, "edgeR") \%>\%
  make_voom(., my_design) \%>\%
  model_limma() \%>\%
  make_contrasts(design_matrix = my_design, Statuspregnant, Statusvirgin) \%>\%
  model_bayes() \%>\%
  tidy.marray.lm() \%>\%
  volcano_plot()
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_analysis.R
\name{make_design_matrix}
\alias{make_design_matrix}
\title{Create a design matrix}
\usage{
make_design_matrix(metadata, vars)
}
\arguments{
\item{metadata}{cleaned metadata for RNAseq data}

\item{vars}{a character vector of variables to include in the model
these should be case-specific and match the variable names in the metadata}
}
\value{
a `tbl` of the design matrix. Column names will be formatted as
`paste0(variablename, samplename)`
}
\description{
`make_design_matrix()` creates a model matrix for your DE analysis.
See \code{\link[stats]{model.matrix}} for further details on what a design
(or model) matrix is. This function prints out the column names so that you
can view what your design matrix variable names are
}
\details{
The order in which you specify your variables will affect which
variable is dummy coded as the reference variable.
}
\examples{
make_design_matrix(metadata, "CellType")

make_design_matrix(metadata, c("CellType", "Status"))

# Results in a different design matrix
make_design_matrix(metadata, c("Status", "CellType"))

# In a pipeline
counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") \%>\%
  mutate(FileName = stringr::str_replace(FileName, "\\\\.", "-"))

my_design <- check_sample_names(counts, c(1,2), meta, FileName) \%>\%
  purrr::pluck("meta") \%>\%
  make_design_matrix(., c("CellType"))
}

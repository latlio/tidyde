% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_analysis.R
\name{check_sample_names}
\alias{check_sample_names}
\title{Functions for simulating RNAseq DE data}
\usage{
check_sample_names(count_df, cols_to_remove, metadata, metadata_var)
}
\arguments{
\item{count_df}{cleaned dataframe of counts, rows should be gene IDs,
columns should be samples, cells should only contain counts}

\item{cols_to_remove}{vector of column numbers that do not correspond to a sample,
necessary to identify for downstream functions}

\item{metadata}{cleaned metadata for RNAseq data}

\item{metadata_var}{column of sample identifier that user expects to match with
count_matrix}
}
\value{
a `list` with the following components:
\item{old_count}{the original count dataframe supplied}
\item{mod_count}{the pure count dataframe (no other columns)}
\item{meta}{sorted metadata (if necessary), otherwise the supplied
metadata is returned with console message output of the quality control check.}
}
\description{
Check sample names
}
\details{
`check_sample_names()` is a simple quality control step that verifies whether
the column names in the count matrix match with a user-defined metadata column.
In order for a match to occur, the value levels of the column names and those of the
user-defined metadata column need to be identical
(e.g. `setdiff(colnames(my_count_matrix), metadata$my_column)` should be 0),
and the order in which the values appear need to be identical
(e.g. `identical(colnames(my_count_matrix), metadata$my_column)` should be TRUE).

The proportion of zeros in the original count data is also printed to the
console. For count data that has a medium to high proportion of zeros,
\code{\link[edgeR]{voomLmFit}} is recommended. Otherwise
\code{\link[edgeR]{voom}} followed by \code{\link[edgeR]{lmFit}} is
recommended.
}
\examples{
counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") \%>\%
mutate(FileName = stringr::str_replace(FileName, "\\\\.", "-"))

check_sample_names(counts, c(1,2), meta, FileName)

}

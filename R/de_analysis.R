#' Check sample names
#'
#' `check_sample_names()` is a simple quality control step that verifies whether
#' the column names in the count matrix match with a user-defined metadata column.
#' In order for a match to occur, the value levels of the column names and those of the
#' user-defined metadata column need to be identical
#' (e.g. `setdiff(colnames(my_count_matrix), metadata$my_column)` should be 0),
#' and the order in which the values appear need to be identical
#' (e.g. `identical(colnames(my_count_matrix), metadata$my_column)` should be TRUE).
#'
#' @param count_df cleaned dataframe of counts, rows should be gene IDs,
#' columns should be samples, cells should only contain counts
#' @param cols_to_remove vector of column numbers that do not correspond to a sample,
#' necessary to identify for downstream functions
#' @param metadata cleaned metadata for RNAseq data
#' @param metadata_var column of sample identifier that user expects to match with
#' count_matrix
#'
#' @return a `list` of `old_count`, `mod_count`, and `meta`
#' where `old_count` is the original count dataframe supplied,
#' `mod_count` is the pure count dataframe (no other columns),
#' and `meta` is the sorted metadata (if necessary), otherwise the supplied
#' metadata is returned with console message output of the quality control check
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#' mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' check_sample_names(counts, EntrezGeneID, c(1,2), meta, FileName)
#'
check_sample_names <- function(count_df,
                               cols_to_remove,
                               metadata,
                               metadata_var) {

  # coerce to dataframe
  if(!is.data.frame(count_df)) {
    message("Your input data was coerced into a data frame.")
    count_df <- as.data.frame(count_df)
  }

  #output count dataframe with purely counts
  count_df_mod <- count_df[,-cols_to_remove]

  if (
    identical(
      setdiff(colnames(count_df)[-cols_to_remove],
              metadata %>% dplyr::select({{metadata_var}}) %>% dplyr::pull()),
              character(0)
    )
  ){
    message("The column names of the count matrix and the unique sample ID values are correctly specified.")

    if(
      identical(
        colnames(count_df)[-cols_to_remove],
        metadata %>% dplyr::select({{metadata_var}}) %>% dplyr::pull()
      )
    ) {
      message("The order is also correct. You can safely proceed with the remaining analysis steps.")
    } else {
      message("The data was ordered incorrectly. The metadata has been reordered.")
      metadata <- metadata[
        match(
          colnames(count_df)[-cols_to_remove],
          metadata %>% dplyr::select({{metadata_var}}) %>% dplyr::pull()
        ),]
    }
  } else {
    stop("Please specify a different variable or check that the values in the metadata are written correctly.")
  }
  list(
    old_count = count_df,
    mod_count = count_df_mod,
    meta = metadata
  )
}

#' Create a design matrix
#'
#' `make_design_matrix()` creates a model matrix for your DE analysis.
#' See \code{\link[stats]{model.matrix}} for further details on what a design
#' (or model) matrix is.
#'
#' @param metadata cleaned metadata for RNAseq data
#' @param vars a character vector of variables to include in the model
#'
#' @details The order in which you specify your variables will affect which
#' variable is dummy coded as the reference variable.
#'
#' @return a `tbl` of the design matrix
#'
#' @export
#'
#' @examples
#'
#' make_design_matrix(metadata, "CellType")
#'
#' make_design_matrix(metadata, c("CellType", "Status"))
#'
#' # Results in a different design matrix
#' make_design_matrix(metadata, c("Status", "CellType"))
#'
#' # In a pipeline
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' my_design <- check_sample_names(counts, EntrezGeneID, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta") %>%
#'   make_design_matrix(., c("CellType"))

make_design_matrix <- function(metadata, vars) {
  formula <- as.formula(paste0("~ 0 + ", paste(vars, collapse = "+")))
  modelr::model_matrix(metadata,
                       formula)
}

#' Filter lowly expressed genes
#'
#' `filter_genes()` is a wrapper function for several filtering methods.
#'
#' @param count_df preprocessed dataframe of pure counts
#' @param id vector of gene IDs
#' @param filter_method Either `edgeR`, `samplenr`, or `cpm`
#' @param min_samples minimum number of samples
#' @param min_cpm minimum cpm
#' @param ... additional arguments to `filterByExpr()`
#'
#' @details I encourage users to exercise caution before using this filter function.
#' Oftentimes, the filtering step should be specific to the sequencing experiment.
#' The `edgeR` option is a wrapper for `edgeR::filterByExpr()`.
#' The `samplenr` option filters out genes across sample whose counts are lower
#' 2*number_of_samples
#' The `cpm` option filters out genes whose rowsums (excluding cells lower
#' than `min_cpm`) are less than number_of_samples/min_samples
#'
#' @return a `list` of a vector of filtered counts called `counts`,
#' a dataframe called `samples` containing the library sizes and the normalization
#' factors, and a dataframe called `genes` contained the gene IDs
#'
#' @export
#'
#' @examples
#'
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' # this step may differ depending on how your data is formatted
#' id <- counts$EntrezGeneId
#'
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR")

filter_genes <- function(count_df,
                         id,
                         filter_method,
                         min_samples = 10,
                         min_cpm = 0.25,
                         ...) {

  filter_by_edgeR <- function(dge, id) {
    counts <- dge$counts
    #rownames is not "tidy" but it may be easier to work with
    rownames(counts) <- id
    genes_to_keep <- edgeR::filterByExpr(dge, ...)
    counts <- counts[genes_to_keep,]
    genes <- rownames(counts)
    dge_filtered <- edgeR::DGEList(counts = counts,
                                   genes = genes)
    dge_filtered
  }

  filter_by_samplenr <- function(dge, id) {
    counts <- dge$counts
    rownames(counts) <- id
    genes_to_keep <- rowSums(counts) >= 2*ncol(counts)
    counts <- counts[genes_to_keep,]
    genes <- rownames(counts)
    dge_filtered <- edgeR::DGEList(counts = counts,
                                   genes = genes)
    dge_filtered
  }

  filter_by_cpm <- function(dge,
                            id,
                            min_samples,
                            min_cpm) {
    counts <- dge$counts
    rownames(counts) <- id
    tmp <- edgeR::cpm(counts) %>% as.data.frame()

    genes_to_keep <- tmp %>%
      tibble::rownames_to_column("id") %>%
      rowwise() %>%
      mutate(cond = {ind <- c_across(2:last_col())
                     sum(ind >= min_cpm) >= ncol(tmp)/min_samples}) %>%
      ungroup() %>%
      tibble::column_to_rownames("id") %>%
      filter(cond)

    counts <- counts[rownames(genes_to_keep),]
    genes <- rownames(counts)
    dge_filtered <- edgeR::DGEList(counts = counts,
                                   genes = genes)
    dge_filtered
  }

  dge <- edgeR::DGEList(count_df)
  dge_filtered <- switch(filter_method,
                         edgeR = filter_by_edgeR(dge, id),
                         samplenr = filter_by_samplenr(dge, id),
                         cpm = filter_by_cpm(dge,
                                             id,
                                             min_samples,
                                             min_cpm))
  dgelist_w_normfactors <- calcNormFactors(dge_filtered)
  dgelist_w_normfactors
}

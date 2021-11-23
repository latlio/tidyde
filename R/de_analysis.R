#' Functions for simulating RNAseq DE data
#'
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
#' @return a `list` with the following components:
#' \item{old_count}{the original count dataframe supplied}
#' \item{mod_count}{the pure count dataframe (no other columns)}
#' \item{meta}{sorted metadata (if necessary), otherwise the supplied
#' metadata is returned with console message output of the quality control check.}
#'
#' @details The proportion of zeros in the original count data is also printed to the
#' console. For count data that has a medium to high proportion of zeros,
#' \code{\link[edgeR]{voomLmFit}} is recommended. Otherwise
#' \code{\link[edgeR]{voom}} followed by \code{\link[edgeR]{lmFit}} is
#' recommended.
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#' mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' check_sample_names(counts, c(1,2), meta, FileName)
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
  if(is.null(cols_to_remove)) {
    count_df_mod <- count_df
  } else {
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
  }

  # cat("The proportion of zeroes in your count data is ",
  #     sum(count_df == 0)/(ncol(count_df) * nrow(count_df)))

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
#' (or model) matrix is. This function prints out the column names so that you
#' can view what your design matrix variable names are
#'
#' @param metadata cleaned metadata for RNAseq data
#' @param vars a character vector of variables to include in the model
#' these should be case-specific and match the variable names in the metadata
#'
#' @details The order in which you specify your variables will affect which
#' variable is dummy coded as the reference variable.
#'
#' @return a `tbl` of the design matrix. Column names will be formatted as
#' `paste0(variablename, samplename)`
#'
#' @export
#'
#' @examples
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
#' my_design <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta") %>%
#'   make_design_matrix(., c("CellType"))

make_design_matrix <- function(metadata, vars) {
  formula <- as.formula(paste0("~ 0 + ", paste(vars, collapse = "+")))
  out <- modelr::model_matrix(metadata,
                              formula)
  print(tibble(variables = colnames(out)))
  out
}

#'
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
#' @return a `list` (`DGEList`) with the following components:
#' \item{counts}{a vector of the filtered counts}
#' \item{samples}{a dataframe containing the library sizes and the normalization
#' factors}
#' \item{genes}{a dataframe containing the gene IDs}
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' # this step may differ depending on how your data is formatted
#' id <- as.character(counts$EntrezGeneID)
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

#' Convert count data to voom
#'
#' `make_voom()` is a wrapper function for `voom()`
#'
#' @param .dge a `list` or a `DGElist` object
#' @param design_matrix a design matrix with rows corresponding to samples
#' and columns to coefficients to be estimated
#' @param .f limma::voom
#' @param ... additional arguments passed to `.f`
#'
#' @details Please refer to \code{\link[limma]{voom}} for more information.
#'
#' @return a `list` with the following components:
#' \item{E}{a numeric matrix of normalized expression values on the log2 scale}
#' \item{weights}{numeric matrix of inverse variance weights}
#' \item{design}{design matrix}
#' \item{lib.size}{numeric vector of total normalized library sizes}
#' \item{genes}{dataframe of gene annotation extracted from `counts`}
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design)

make_voom <- function(.dge, design_matrix, .f = limma::voom, ...) {
  .args <- rlang::enexprs(...)

  rlang::eval_tidy(rlang::expr(.f(counts = .dge,
                                  design = design_matrix,
                                  !!! .args)))
}

#' Model limma
#'
#' `model_limma()` is a wrapper function for `lmFit`. Fits a linear model for
#' each gene.
#'
#' @param voom a voom object
#' @param .f limma::lmFit
#' @param ... additional arguments to .f
#'
#' @details Please refer to \code{\link[limma]{lmFit}} for more information.
#'
#' @return a `list` (`MArrayLM`) object containing the results of the fit
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma()

model_limma <- function(.data, .f = limma::lmFit, ...) {
  .args <- rlang::enexprs(...)

  rlang::eval_tidy(rlang::expr(.f(object = .data,
                                  !!! .args)))
}

#' Make contrasts
#'
#' `make_contrasts()` computes estimated coefficients and standard errors
#' for a given set of contrasts
#'
#' @param .fit an MArrayLM object or list object. Must contain components
#' `coefficients` and `stdev.unscaled`
#' @param design_matrix a design matrix with rows corresponding to samples
#' and columns to coefficients to be estimated
#' @param .condition1 the control
#' @param .condition2 the experimental
#' @param ... the names of the variables which you like to compare
#'
#' @details Please refer to \code{\link[limma]{contrasts.fit}} for more information.
#'
#' @return a `list` object of the same class as .fit
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(Statuspregnant, Statusvirgin)
make_contrasts <- function(.fit, design_matrix, .condition1, .condition2) {
  # check types, if input is symbol
  if(is.symbol(.condition1)) {
    condition1 <- rlang::enexpr(.condition1)
  } else {
    #if string
    condition1 <- .condition1
  }

  if(is.symbol(.condition2)) {
    condition2 <- rlang::enexpr(.condition2)
  } else {
    #if string
    condition2 <- .condition2
  }

  my_contrast <- paste(c(condition2, condition1), collapse = " - ")
  contrast_matrix <- limma::makeContrasts(
    contrasts = c(my_contrast),
    levels = colnames(design_matrix)
  )
  contrasts.fit(.fit, contrast_matrix)
}

#' Model differential expression
#'
#' `model_bayes()` performs an empirical Bayes fit
#'
#' @param .fit an MArrayLM object produced by `model_limma()` or
#' `make_contrasts()`
#' @param .f limma::eBayes
#' @param ... additional arguments to .f
#'
#' @details Please refer to \code{\link[limma]{eBayes}} for more information.
#'
#' @return a `list` (`MArrayLM`) object containing the results of the fit
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(Statuspregnant, Statusvirgin) %>%
#'   model_bayes()
model_bayes <- function(.fit, .f = limma::eBayes, ...) {
  .args <- rlang::enexprs(...)

  rlang::eval_tidy(rlang::expr(.f(fit = .fit,
                                  !!! .args)))
}

#' A function that runs a sample limma-voom DE workflow
#'
#' This wrapper function uses all the functions within this function and strings them
#' together in a sample workflow. Currently only `limma` is supported.
#'
#' @param raw_counts cleaned dataframe of counts, rows should be gene IDs,
#' columns should be samples, cells should only contain counts
#' @param cols_to_remove vector of column numbers that do not correspond to a sample,
#' necessary to identify for downstream functions
#' @param metadata cleaned metadata for RNAseq data
#' @param .metadata_var column of sample identifier that user expects to match with
#' count_matrix
#' @param .design_var a character vector of variables to include in the model
#' these should be case-specific and match the variable names in the metadata
#' @param contrastcond1 control var
#' @param contrastcond2 experimentl var
#' @param .filter_method Please refer to \code{\link[tidyde]{filter_genes}} for more information.
#'
#' @return a `tbl`
#'
#' @export
#'
#' @examples
#'
#' cond_df <- tibble(cond1 = c("a1", "a2", "a3"), cond2 = c("b1", "b2", "b3"))
#'
#' map2(cond_df$cond1, cond_df$cond2, ~run_sample_limma_de(
#' my_counts,
#' my_metadata,
#' my_variable,
#' my_design_variable,
#' .x,
#' .y))
run_sample_limma_de <- function(raw_counts,
                                metadata,
                                .metadata_var,
                                .design_vars,
                                contrastcond1,
                                contrastcond2,
                                .cols_to_remove = NULL,
                                .filter_method = "edgeR"
) {
  the_design <- check_sample_names(raw_counts,
                                   cols_to_remove = .cols_to_remove,
                                   metadata,
                                   .metadata_var) %>%
    purrr::pluck("meta") %>%
    make_design_matrix(., .design_vars)

  the_id <- as.character(row.names(raw_counts))

  the_de_res <- check_sample_names(raw_counts,
                                   cols_to_remove = .cols_to_remove,
                                   metadata,
                                   .metadata_var) %>%
    purrr::pluck("mod_count") %>%
    filter_genes(., the_id, .filter_method) %>%
    make_voom(., the_design) %>%
    model_limma() %>%
    make_contrasts(design_matrix = the_design,
                   contrastcond1,
                   contrastcond2) %>%
    model_bayes()

  tidy.marray.lm(the_de_res)
}

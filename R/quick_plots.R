#' Quickly plot a boxplot of counts
#'
#' `counts_boxplot()` generates a boxplot of counts given count data.
#'
#' @param count_df cleaned dataframe of counts, rows should be gene IDs,
#' columns should be samples, cells should only contain counts
#' @param metadata cleaned metadata for RNAseq data
#' @param sample_var variable of sample ids
#' @param facet_var variable to facet plots by
#' @param .y_intercept default 2000000, setting to NULL removes horizontal line
#'
#' @return a boxplot of RNAseq counts (ggplot object)
#'
#' @details Like any other ggplot object, you can customize the theme of the plot.
#' Note that this function generates a fairly generic boxplot, and thus is aimed
#' for quick exploratory purposes (much like the intention behind `qplot()`).
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#' mycount <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count")
#' mymeta <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta")
#'
#' # View total raw counts
#'
#' counts_boxplot(mycount, mymeta, SampleName)
#'
#' # View filtered counts and also facet by variable
#' id <- as.character(counts$EntrezGeneID)
#' mycount %>%
#'   filter_genes(., id, "edgeR") %>%
#'   purrr::pluck("counts") %>%
#'   counts_boxplot(., mymeta, SampleName, CellType)
counts_boxplot <- function(count_df, metadata, sample_var, facet_var = NULL,
                           .y_intercept = 2000000) {

  data_to_plot <- tibble(
    totalcounts = as.vector(colSums(count_df))
  ) %>% bind_cols(metadata)

  if(!rlang::quo_is_null(enquo(facet_var))) {
    ggplot(data_to_plot,
           aes(x = {{sample_var}},
               y = totalcounts,
               fill = {{sample_var}})) +
      geom_boxplot() +
      facet_wrap(vars({{facet_var}}), nrow = 1) +
      theme_bw() +
      ylab("Total counts") + xlab("") +
      geom_jitter(position = position_jitter(0.2),
                  size = 0.8) +
      geom_hline(yintercept = .y_intercept,
                 linetype = "dashed",
                 color = "red")  +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = 1))
  } else {
    ggplot(data_to_plot,
           aes(x = {{sample_var}},
               y = totalcounts,
               fill = {{sample_var}})) +
      geom_boxplot() +
      theme_bw() +
      ylab("Total counts") + xlab("") +
      geom_jitter(position = position_jitter(0.2),
                  size = 0.8) +
      geom_hline(yintercept = .y_intercept,
                 linetype = "dashed",
                 color = "red")  +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = 1))
  }
}
#'
#' Quickly plot a PCA plot
#'
#' `pca_plot()`.
#'
#' Quickly plot a volcano plot
#'
#' `volcano_plot()`

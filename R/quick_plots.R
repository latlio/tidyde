#' Quickly plot a boxplot of counts
#'
#' `counts_boxplot()` generates a boxplot of counts given count data.
#'
#' @param count_df cleaned dataframe of counts, rows should be gene IDs,
#' columns should be samples, cells should only contain counts. Also accepts a
#' `DGEList` object
#' @param metadata cleaned metadata for RNAseq data
#' @param sample_var variable of sample ids
#' @param facet_var variable to facet plots by
#' @param .y_intercept default 2000000, setting to NULL removes horizontal line
#'
#' @return a boxplot of RNAseq counts (ggplot object)
#'
#' @details Like any other ggplot object, you can customize the theme of the plot.
#' Note that this function generates a fairly generic boxplot, and thus is intended
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
#' counts_boxplot(mycount, mymeta, SampleName)
#'
#' # View filtered counts and also facet by variable
#' id <- as.character(counts$EntrezGeneID)
#' mycount %>%
#'   filter_genes(., id, "edgeR") %>%
#'   counts_boxplot(., mymeta, SampleName, CellType)
#'
#' # As an EDA step within a pipeline of functions
#' my_design <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta") %>%
#'   make_design_matrix(., c("Status"))
#'
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %T>%
#'   {print(counts_boxplot(., check_sample_names(counts, c(1,2), meta, FileName) %>%
#'     purrr::pluck("meta"), SampleName))} %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(design_matrix = my_design, Statuspregnant, Statusvirgin) %>%
#'   model_bayes()
counts_boxplot <- function(count_df, metadata, sample_var, facet_var = NULL,
                           .y_intercept = 2000000) {

  if(inherits(count_df, "DGEList")) {
    count_df <- count_df[["counts"]]
  }

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
#' `pca_plot()` plots the first two PC in a scatterplot.
#'
#' @param voom_data, voom object
#' @param metadata cleaned metadata for RNAseq data
#' @param ... arguments to aes
#'
#' @return a scatterplot of the first two PC (ggplot object)
#'
#' @details Like any other ggplot object, you can customize the theme of the plot.
#' Note that this function generates a fairly generic PCA plot, and thus is intended
#' for quick exploratory purposes (much like the intention behind `qplot()`).
#'
#' @export
#'
#' @examples
#' # As an EDA step within a pipeline of functions
#'
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' my_design <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta") %>%
#'   make_design_matrix(., c("Status"))
#'
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %T>%
#'   {print(pca_plot(., check_sample_names(counts, c(1,2), meta, FileName) %>%
#'     purrr::pluck("meta")))} %>%
#'   model_limma() %>%
#'   make_contrasts(design_matrix = my_design, Statuspregnant, Statusvirgin) %>%
#'   model_bayes()
pca_plot <- function(voom_data, metadata, ...) {

  calculate_pca_var_explained <- function(prcomp_obj,
                                          index) {
    ## Function Purpose: get the variance explained by PCA component i

    ## Inputs:
    #' prcomp_obj, result from prcomp()
    #' index, which PCA you want
    eigs <- prcomp_obj$sdev^2
    eigs[index] / sum(eigs)
  }

  pca <- prcomp(t(voom_data$E), scale = TRUE, tol = 0)
  pca_data <- as_tibble(pca$x) %>%
    bind_cols(metadata %>%
                mutate(across(everything(), as.factor)))

  ggplot(pca_data, aes_string(...)) +
    aes(x = PC1, y = PC2) +
    geom_point(size = 3) +
    theme_bw() +
    labs(x = glue::glue("PCA 1 (",
                        {round(calculate_pca_var_explained(pca, 1) * 100,
                               1)},
                        "%)"),
         y = glue::glue("PCA 2 (",
                        {round(calculate_pca_var_explained(pca, 2) * 100,
                               1)},
                        "%)"))
  }
#' Quickly plot a volcano plot
#'
#' `volcano_plot()` plots a volcano plot
#'
#' @param fit tidied dataframe of fit results
#' @param fc_cutoff cutoff to determine what is considered differentially expressed
#' @param pval_cutoff plots horizontal line at p-value
#'
#' @return a `list` containing a volcano plot (ggplot object) and the data used to
#' generate the volcano plot
#'
#' @details Like any other ggplot object, you can customize the theme of the plot.
#' Note that this function generates a fairly generic volcano plot, and thus is
#' intended for quick exploratory purposes (much like the intention behind `qplot()`).
#'
#' @export
#'
#' @examples
#' # As an EDA step within a pipeline of functions
#'
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' my_design <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta") %>%
#'   make_design_matrix(., c("Status"))
#'
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(design_matrix = my_design, Statuspregnant, Statusvirgin) %>%
#'   model_bayes() %>%
#'   tidy.marray.lm() %>%
#'   volcano_plot()

volcano_plot <- function(fit, fc_cutoff = 2, pval_cutoff = 0.05) {

  data_to_plot <- fit %>%
    dplyr::mutate(DE = case_when(
      estimate > fc_cutoff & p.value < pval_cutoff ~ "Up",
      estimate < -1 * fc_cutoff & p.value < pval_cutoff ~ "Down",
      TRUE ~ "Not DE"),
      color_groups = case_when(
        p.value >= pval_cutoff ~ "NS",
        p.value < pval_cutoff & abs(estimate) < fc_cutoff ~ "Significant",
        p.value >= pval_cutoff & abs(estimate) >= fc_cutoff ~ "High FC",
        p.value < pval_cutoff & abs(estimate) >= fc_cutoff ~ "Significant and high FC"),
      alpha = case_when(
        DE == "Up" ~ 1,
        DE == "Down" ~ 1,
        TRUE ~ 0.3
      ))

  outp <- ggplot(data = data_to_plot,
                 aes(x = estimate,
                     y = -log10(p.value),
                     col = color_groups,
                     label = gene,
                     alpha = alpha)) +
    geom_point() +
    geom_label_repel(data = subset(data_to_plot,
                                   color_groups == "Significant and high FC"),
                     aes(label = gene),
                     col = "black",
                     alpha = 1,
                     show.legend = FALSE,
                     min.segment.length = unit(0, 'lines')) +
    geom_segment(
      x = 0.5,
      y = max(-log10(data_to_plot$p.value), na.rm = TRUE) + 4,
      xend = max(data_to_plot$estimate, na.rm = TRUE) + 1,
      yend = max(-log10(data_to_plot$p.value), na.rm = TRUE) + 4,
      lineend = "square", # See available arrow types in example above
      linejoin = "round",
      size = 1,
      arrow = arrow(length = unit(0.3, "inches")),
      colour = "#31938A") +
    annotate(geom = "label",
             x = max(data_to_plot$estimate, na.rm = TRUE) + 1,
             y = max(-log10(data_to_plot$p.value), na.rm = TRUE) + 4.5,
             label = "Treatment",
             col = "#31938A") +
    geom_segment(
      x = -0.5,
      y = max(-log10(data_to_plot$p.value), na.rm = TRUE) + 4,
      xend = min(data_to_plot$estimate, na.rm = TRUE) - 1,
      yend = max(-log10(data_to_plot$p.value), na.rm = TRUE) + 4,
      lineend = "square", # See available arrow types in example above
      linejoin = "round",
      size = 1,
      arrow = arrow(length = unit(0.3, "inches")),
      colour = "#31938A") +
    annotate(geom = "label",
             x = min(data_to_plot$estimate, na.rm=TRUE) - 1,
             y = max(-log10(data_to_plot$p.value), na.rm = TRUE) + 4.5,
             label = "Control",
             col = "#31938A") +
    geom_vline(xintercept = c(-1 * fc_cutoff, fc_cutoff),
               lty = 2, alpha = 0.2) +
    geom_hline(yintercept = -log10(pval_cutoff),
               lty = 2, alpha = 0.2) +
    xlim(min(data_to_plot$estimate, na.rm = TRUE) - 1.5,
         max(data_to_plot$estimate, na.rm = TRUE) + 1.5) +
    ylim(0, max(-log10(data_to_plot$p.value), na.rm = TRUE) + 5) +
    labs(x = bquote(~-Log[2]~ "fold change"),
         y = bquote(~-Log[10]~adjusted~italic(P))) +

    scale_color_manual(values = c("black", "#4E91FF", "#008000", "#FE6C61")) +
    scale_alpha(guide = "none")

    list(data = data_to_plot,
         volcano = outp)
}

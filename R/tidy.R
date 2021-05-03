#' Tidying methods for MArrayLM objects
#'
#' These methods tidy the results of differential expression analysis objects.
#' Currently only `limma` is supported.
#'
#' @param x a differential expression fit object. Currently supports `MArrayLM`
#' from the `limma` package.
#' @param conf.int logical. Include confidence intervals?
#' @param exponentiate logical. Should the estimates and (if `conf.int` = `TRUE`)
#' confidence intervals be exponentiated?
#' @param ... additional arguments
#'
#' @return a `tbl`
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' out <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(Statuspregnant, Statusvirgin) %>%
#'   model_bayes()
#'
#' out %>% tidy.marray.lm()
#'
#' @rdname tidiers
tidy.marray.lm <- function(x, conf.int = TRUE, exponentiate = FALSE, ...) {
  if (!inherits(x, "MArrayLM")) stop("`x` must be of class `MArrayLM`")

  x <- purrr::map(x, ~as.vector(.x))
  margin_error <- sqrt(x$s2.post)*
    x$stdev.unscaled*qt(0.975, df = x$df.total)
  results <- tibble(gene = unlist(x$genes),
                    estimate = x$coefficients,
                    std.error = x$s2.post,
                    statistic = x$t,
                    p.value = x$p.value,
                    conf.low = x$coefficients - margin_error,
                    conf.high = x$coefficients + margin_error,
                    stringsAsFactors = FALSE)

  if (exponentiate) {
    results$estimate <- exp(results$estimate)
    results$conf.low <- exp(results$conf.low)
    results$conf.high <- exp(results$conf.high)
  }

  if (!conf.int) {
    results <- dplyr::select(results, -conf.low, -conf.high)
  }

  results
}

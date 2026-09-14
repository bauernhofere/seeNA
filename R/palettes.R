#' Default copy-number state colors
#'
#' The named palette used by [plot_ichor_profile()] for the four display
#' states returned by [ichor_call_state()].
#' @return A named character vector.
#' @examples
#' ichor_state_colors()
#' @export
ichor_state_colors <- function() {
  c("Deep loss" = "#A6E4C0", "Loss" = "#1E9E5E",
    "Neutral" = "#2E6CB0", "Gain" = "#CB2B2B")
}

.default_sample_colors <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

.theme_ichor <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "top"
    )
}

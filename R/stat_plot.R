#' Calculate Internal Flow Statistics (Optimized)
#'
#' Computes summary statistics for internal relative flows from a list of matrices.
#' This optimized version uses efficient base R functions and is well-documented.
#'
#' @param matrix_list A list of square numeric matrices. The list must be named
#'   with the round numbers, which correspond to the number of clusters.
#' @return A data frame with columns: `round`, `mean`, `min`, `median`, `max`.
#' @export
flowbca_stat <- function(matrix_list) {
  # sapply is an efficient way to iterate over a list and collect results.
  stats_matrix <- sapply(matrix_list, function(mat) {
    row_sums <- rowSums(mat)
    # Handle cases where row sum is 0 to prevent division by zero (NaN/Inf).
    internal_flows <- ifelse(row_sums == 0, 0, diag(mat) / row_sums)

    # Return a named vector of summary statistics for each matrix.
    c(
      mean = mean(internal_flows, na.rm = TRUE),
      min = min(internal_flows, na.rm = TRUE),
      median = median(internal_flows, na.rm = TRUE),
      max = max(internal_flows, na.rm = TRUE)
    )
  })

  # The result from sapply is a matrix with stats in rows and rounds in columns.
  # Transpose it to have rounds in rows, then convert to a data frame.
  stats_df <- as.data.frame(t(stats_matrix))

  # Add round numbers from the list names as a proper column.
  stats_df$round <- as.numeric(rownames(stats_df))

  # Ensure the column order is logical.
  stats_df <- stats_df[, c("round", "mean", "min", "median", "max")]
  rownames(stats_df) <- NULL # Reset row names for a clean data frame.

  return(stats_df)
}

#' Plot Flow Statistics Over Clustering Rounds (Optimized)
#'
#' Visualizes flow statistics with a fixed y-axis (0-1) and dynamic plot styles.
#' This optimized version is more robust, readable, and handles edge cases.
#'
#' @param stat_data A data frame produced by `flowbca_stat`.
#' @param upper_bound Optional: the maximum number of rounds (clusters) to display.
#' @return A plot is drawn on the current graphics device. Returns `invisible(NULL)`.
#' @export
flowbca_plot <- function(stat_data, upper_bound = NULL) {
  # --- 1. Input Validation ---
  required_cols <- c("round", "mean", "min", "median", "max")
  if (!all(required_cols %in% names(stat_data))) {
    stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(stat_data$round)) {
    stat_data$round <- as.numeric(as.character(stat_data$round))
  }

  # --- 2. Data Filtering ---
  plot_df <- stat_data
  if (!is.null(upper_bound) && is.numeric(upper_bound)) {
    plot_df <- plot_df[plot_df$round <= upper_bound, ]
  }
  if (nrow(plot_df) == 0) {
    warning("No data to plot after applying the filter.")
    return(invisible(NULL))
  }
  plot_df <- plot_df[order(plot_df$round, decreasing = TRUE), ]

  # --- 3. Dynamic Plot Style & Edge Case Handling ---
  num_points <- nrow(plot_df)
  plot_type <- if (num_points <= 20) "b" else "l"
  plot_pch <- if (num_points <= 20) 19 else NA
  x_values <- plot_df$round
  xlim_range <- if (num_points > 1) {
    c(max(x_values, na.rm = TRUE), min(x_values, na.rm = TRUE))
  } else {
    x_values + c(0.5, -0.5) # Add padding for a single point plot
  }

  # --- 4. Plotting ---
  y_values <- plot_df[, c("mean", "min", "median", "max")]
  plot_colors <- c("blue", "red", "green", "purple")

  matplot(x = x_values, y = y_values, type = plot_type, pch = plot_pch,
          lty = 1, lwd = 2, col = plot_colors, xlim = xlim_range, ylim = c(0, 1),
          xaxt = "n", # Suppress default x-axis
          xlab = "Round (Number of Clusters + 1)", ylab = "Internal Relative Flow",
          main = "Flow Statistics per Clustering Round")
  grid() # Add a grid for better readability

  # Add custom integer-only x-axis
  axis_ticks <- pretty(x_values)
  # Filter to keep only integer ticks
  axis_ticks <- axis_ticks[axis_ticks == floor(axis_ticks)]
  axis(1, at = axis_ticks, labels = axis_ticks)

  # --- 5. Add Text Labels (if few points) ---
  if (num_points <= 20) {
    stat_names <- c("mean", "min", "median", "max")
    # Use mapply for a more idiomatic R way to loop over multiple vectors
    mapply(function(y_col, color) {
      y_vals <- plot_df[[y_col]]
      text(x = x_values, y = y_vals, labels = round(y_vals, 3),
           col = color, pos = 3, cex = 0.75)
    }, stat_names, plot_colors)
  }

  # --- 6. Add Legend ---
  legend("topleft", legend = names(y_values), col = plot_colors, lwd = 2, bty = "n")

  return(invisible(NULL))
}
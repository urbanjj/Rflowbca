#' Find Hierarchy Path using Base R (Internal Helper)
#'
#' Traverses a parent-child map to find the full hierarchy path for a single node.
#' This is an internal helper function implemented using only base R.
#'
#' @param start_node The node to start the traversal from.
#' @param cluster_map A named vector where names are children and values are parents.
#' @return A string representing the full path.
find_hierarchy <- function(start_node, cluster_map) {
  path <- c(start_node)
  current_node <- start_node
  visited <- c(start_node)

  while (current_node %in% names(cluster_map)) {
    parent_node <- cluster_map[[current_node]]

    # Stop if the parent is NA or empty
    if (is.na(parent_node)) {
      break
    }

    # Check for circular references
    if (parent_node %in% visited) {
      warning(paste("Circular reference detected:", paste(c(visited, parent_node), collapse=" -> ")))
      path <- c(path, "...") # Indicate a circular path
      break
    }

    path <- c(path, parent_node)
    visited <- c(visited, parent_node)
    current_node <- parent_node
  }

  return(paste(rev(path), collapse = "/"))
}

#' Build and Add Hierarchy Path Column using Base R
#'
#' This function constructs hierarchy paths from a data frame of parent-child
#' relationships using only base R functions. It is a replacement for the
#' original `find_hierarchy_db`.
#'
#' @param data A data frame, typically the `unit_set` from a `flowbca` result.
#' @param child_col A string, the name of the column with child units (defaults to "sourceunit").
#' @param parent_col A string, the name of the column with parent units (defaults to "destinationunit").
#' @return The input data frame with a new `hierarchy` column.
#' @importFrom stats setNames
#' @export
build_hierarchy <- function(data, child_col = "sourceunit", parent_col = "destinationunit") {
  # Ensure input is a data frame
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  # Create the parent-child map from the data frame
  valid_links <- data[!is.na(data[[parent_col]]), ]
  cluster_map <- setNames(valid_links[[parent_col]], valid_links[[child_col]])
  clusterid <- unique(data$clusterid)

  # Use sapply to apply the find_hierarchy function to each child unit
  data$hierarchy <- sapply(data[[child_col]], find_hierarchy, cluster_map = cluster_map)
  data$h_level <- sapply(data$hierarchy, \(x) nchar(x)-nchar(gsub('/','',x)))+1
  data$h_parent <- sapply(data$sourceunit, \(x) ifelse(x %in% clusterid,1,0))
  return(data)
}

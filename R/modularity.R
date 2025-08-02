#' Calculate Modularity for each step of flowbca
#'
#' This function calculates the modularity for each step of the
#' Flow based-Clustering Algorithm (flowbca) clustering process.
#'
#' @param unit_set A data frame from the result of flowbca function.
#'        It must contain 'sourceunit', 'destinationunit', and 'round' columns.
#' @param F_matrix_history A list of flow matrices from the flowbca process.
#'        The first matrix in the list (the original OD matrix) is used for the calculation.
#' @return A data frame with two columns: 'round' and 'modularity',
#'         showing the modularity value for each clustering step.
#' @export
#' @examples
#' # Assuming 'bca_result' is an object from a flowbca function call
#' # modularity_result <- flowbca_modularity(bca_result$unit_set, bca_result$F_matrix_history)
#' # plot(rev(modularity_result$round$round), modularity_result$round$modularity)
flowbca_modularity <- function(unit_set, F_matrix_history){
  # --- Input Validation ---
  if (!is.data.frame(unit_set) || !all(c("sourceunit", "destinationunit", "round") %in% names(unit_set))) {
    stop("`unit_set` must be a data frame with 'sourceunit', 'destinationunit', and 'round' columns.")
  }
  if (!is.list(F_matrix_history) || length(F_matrix_history) == 0) {
    stop("`F_matrix_history` must be a non-empty list of matrices.")
  }

  # --- Main Logic ---
  F_matrix_1 <- F_matrix_history[[1]]
  cluster_history <- cluster_set(unit_set[match(colnames(F_matrix_1),unit_set$sourceunit),])
  M_matrix <- CalcModMatrix(F_matrix_1, is_directed = TRUE, modularity_resolution = 1)

  v_modularity <- unlist(lapply(cluster_history, function(x) {
    calculate_modularity(F_matrix_1, M_matrix, x$clusterid)
  }))

  modularity_df <- data.frame(
    round = names(v_modularity),
    modularity = v_modularity,
    row.names = NULL
  )
  return(modularity_df)
}

#' Create a set of cluster assignments for each merge event.
#' @noRd
cluster_set <- function(unit_set){
  unit_set <- as.data.frame(unit_set)
  round_set <- unit_set[!is.na(unit_set$round), ]
  round_set <- round_set[order(round_set$round, decreasing = TRUE), ]

  cluster_history <- list()
  cluster_history[[1]] <- data.frame(sourceunit = unit_set$sourceunit,
                                     clusterid = unit_set$sourceunit,
                                     stringsAsFactors = FALSE)
  rd <- round_set$round
  # Iterate through each merge event by its row index
  for(i in 1:nrow(round_set)) {
    merge_info <- round_set[i, ] # Ensures only one merge is processed
    previous_clusters <- cluster_history[[i]]
    new_clusters <- previous_clusters
    new_clusters$clusterid[new_clusters$clusterid == merge_info$sourceunit] <- merge_info$destinationunit
    cluster_history[[i+1]] <- new_clusters
  }

  max_rd <- if(length(rd) > 0) max(rd) else 0
  names(cluster_history) <- c(as.character(max_rd + 1), as.character(rd))
  return(cluster_history)
}


#' Calculate the modularity matrix.
#' @noRd
CalcModMatrix <- function(F_matrix_1, is_directed = TRUE, modularity_resolution = 1.0){
  if(is_directed == TRUE) {
    mat <- as.matrix(F_matrix_1)
  } else {
    mat <- as.matrix(F_matrix_1) + t(as.matrix(F_matrix_1))
  }
  total_weight <- sum(mat)
  if (total_weight == 0) {
    return(matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
  }

  out_degree <- rowSums(mat)
  in_degree <- colSums(mat)
  E_degree <- outer(out_degree, in_degree) / total_weight
  m_modularity_matrix <- mat - (modularity_resolution * E_degree)
  return(m_modularity_matrix)
}

#' Internal function to calculate modularity for a given community structure.
#' @noRd
calculate_modularity <- function(F_matrix_1, M_matrix, communities, is_directed = TRUE){
  if(is_directed == TRUE) {
    mat <- as.matrix(F_matrix_1)
  } else {
    mat <- as.matrix(F_matrix_1) + t(as.matrix(F_matrix_1))
  }
  total_weight <- sum(mat)
  same_community_matrix <- outer(communities, communities, "==")
  modularity <- sum(M_matrix[same_community_matrix]) / total_weight
  return(modularity)
}
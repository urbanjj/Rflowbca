#' A Flow-based Cluster Algorithm in R
#'
#' This function is an R translation of the flowbca.ado Stata package
#' by Jordy Meekes (2018). It implements a hierarchical clustering algorithm
#' for flow data, including the detailed tie-breaking rules from the original
#' Mata implementation.
#'
#' @param data A data frame where the first column contains source unit identifiers
#'   and subsequent columns represent the flows to destination units. The destination
#'   columns must be in the same order as the source unit rows.
#' @param q A numeric flow threshold (can be relative or absolute) that will
#'   adjust the default stopping criterion. The algorithm stops if the maximum
#'   flow for merging is below this value. Defaults to 0.
#' @param k An integer specifying the desired number of distinct clusters. The
#'   algorithm stops when this number of clusters is reached. Defaults to 1.
#' @param la A numeric value for the average of the internal relative flows of all
#'   clusters. The algorithm terminates if the calculated average (La) exceeds
#'   this value. Defaults to 1.1 (effectively off).
#' @param lw A numeric value for the weighted average of the internal relative flows
#'   of all clusters. The algorithm terminates if the calculated weighted average
#'   (Lw) exceeds this value. Defaults to 1.1 (effectively off).
#' @param lm A numeric value for the minimum internal relative flow. The algorithm
#'   terminates if the calculated minimum (Lm) exceeds this value. Defaults to
#'   1.1 (effectively off).
#' @param opt_f An integer specifying the optimization function:
#'   1: Directed relative flows (default)
#'   2: Undirected relative flows
#'   3: Directed absolute flows
#'   4: Undirected absolute flows
#' 
#' @param save_k Return list of matrices
#'
#' @return A list containing two data frames:
#'   - `unit_set`: Details of the final cluster assignment for each original unit.
#'   - `cluster_set`: Statistics for the final clusters.
#' 
#' @importFrom stats weighted.mean
#' @export
flowbca <- function(data, q = 0, k = 1, opt_f = 1, la = 1.1, lw = 1.1, lm = 1.1, save_k = FALSE) {

  # --- 1. Initial Setup ---
  source_units <- data[[1]]
  destination_units <- colnames(data)[-1]
  
  if (length(source_units) != length(unique(source_units))) {
    stop("Source unit IDs must be unique names.")
  }

  if (length(destination_units) != length(unique(destination_units))) {
    stop("Destination unit IDs must be unique names.")
  }

  if (!identical(as.character(source_units), destination_units)) {
    stop("Source unit IDs must be identical to and in the same order as destination column names.")
  }

  F_matrix <- as.matrix(data[, -1])
  rownames(F_matrix) <- colnames(F_matrix) <- source_units
  
  if (nrow(F_matrix) != ncol(F_matrix)) {
    stop("The number of source units must equal the number of destination units.")
  }
  if (k < 1 || k > nrow(F_matrix)) {
    stop(paste("k must be between 1 and", nrow(F_matrix)))
  }

  merge_history <- list()
  F_matrix_history <- list()
  
  # --- 2. Main Clustering Loop ---
  while (nrow(F_matrix) > k) {
    
    K <- nrow(F_matrix)
    current_ids <- rownames(F_matrix)

    # --- 2a. Prepare Search Matrix (G_matrix) and Flow Matrix (F_prime) ---
    search_matrix <- NULL
    F_prime <- F_matrix # Base matrix for tie-breaking absolute flows
    
    if (opt_f %in% c(1, 2)) { # Relative flows
      row_sums <- rowSums(F_matrix)
      row_sums[row_sums == 0] <- 1 
      G_matrix <- F_matrix / row_sums
      
      if (opt_f == 1) { # Directed relative
        search_matrix <- G_matrix
      } else { # Undirected relative
        search_matrix <- G_matrix + t(G_matrix)
      }
      F_prime <- G_matrix # Use relative flows for tie-breaking
    } else { # Absolute flows
      if (opt_f == 3) { # Directed absolute
        search_matrix <- F_matrix
      } else { # Undirected absolute
        search_matrix <- F_matrix + t(F_matrix)
      }
    }
    
    diag(search_matrix) <- -Inf 
    
    # --- 2b. Identify Units to Merge (r and s) ---
    max_flow <- max(search_matrix, na.rm = TRUE)
    
    if (max_flow < q) {
      message("Stopping: Maximum flow is below threshold q.")
      break
    }
    
    indices <- which(search_matrix == max_flow, arr.ind = TRUE)
    
    # --- 2c. Disambiguation for non-unique (r,s) pairs (Mata caveats) ---
    if (nrow(indices) > 1) {
      cand_r <- sort(unique(indices[, "row"]))
      cand_s <- sort(unique(indices[, "col"]))
      
      # Use directed flows for tie-breaking in undirected cases
      tie_break_matrix <- if (opt_f %in% c(2, 4)) F_matrix else F_prime
      diag(tie_break_matrix) <- -Inf

      # Caveat 1 & 2: One dimension is unique, the other is not.
      if (length(cand_r) > 1 && length(cand_s) == 1) {
        # Multiple r, one s. Find r with max flow from other r's.
        sub_matrix <- tie_break_matrix[cand_r, cand_r]
        if(sum(sub_matrix) > 0) {
          col_sums <- colSums(sub_matrix)
          r_idx <- cand_r[which.max(col_sums)]
        } else {
          r_idx <- cand_r[1] # Fallback
        }
        s_idx <- cand_s[1]
      } else if (length(cand_r) == 1 && length(cand_s) > 1) {
        # One r, multiple s. Find s with max flow from other s's.
        sub_matrix <- tie_break_matrix[cand_s, cand_s]
        if(sum(sub_matrix) > 0) {
          col_sums <- colSums(sub_matrix)
          s_idx <- cand_s[which.max(col_sums)]
        } else {
          s_idx <- cand_s[1] # Fallback
        }
        r_idx <- cand_r[1]
      } else { # Caveat 3 & 4: Multiple r and multiple s, or complex ties
        # First, try to select the best 's' based on flows between candidate s's
        s_sub_matrix <- tie_break_matrix[cand_s, cand_s]
        if (sum(s_sub_matrix) > 0) {
          s_col_sums <- colSums(s_sub_matrix)
          s_idx <- cand_s[which.max(s_col_sums)]
        } else {
          # Fallback: if no internal flows, just pick the first one
          s_idx <- cand_s[1]
        }
        
        # Now, given the chosen 's', find the best 'r' that flows to it
        possible_r_for_s <- indices[indices[, "col"] == s_idx, "row"]
        if (length(possible_r_for_s) == 1) {
          r_idx <- possible_r_for_s
        } else { # Still multiple r's for our chosen s, apply caveat 1 logic
          r_sub_matrix <- tie_break_matrix[possible_r_for_s, possible_r_for_s]
          if(sum(r_sub_matrix) > 0) {
            r_col_sums <- colSums(r_sub_matrix)
            r_idx <- possible_r_for_s[which.max(r_col_sums)]
          } else {
            r_idx <- possible_r_for_s[1] # Fallback
          }
        }
      }
    } else {
      r_idx <- indices[1, "row"]
      s_idx <- indices[1, "col"]
    }

    # --- 2d. Check Other Stopping Conditions (la, lw, lm) ---
    if (la < 1.1 || lw < 1.1 || lm < 1.1) {
      row_flows_stop <- rowSums(F_matrix)
      internal_flows <- diag(F_matrix) / row_flows_stop
      internal_flows[is.nan(internal_flows) | is.infinite(internal_flows)] <- 0

      La <- mean(internal_flows, na.rm = TRUE)
      Lw <- weighted.mean(internal_flows, row_flows_stop, na.rm = TRUE)
      Lm <- min(internal_flows, na.rm = TRUE)
      
      if (la <= La || lw <= Lw || lm <= Lm) {
        message("Stopping: Condition (la, lw, or lm) met.")
        break
      }
    }
    
    # --- 2e. Aggregate Clusters using Transformation Matrix C ---
    r_id <- current_ids[r_idx]
    s_id <- current_ids[s_idx]
    
    merge_history[[length(merge_history) + 1]] <- list(
      round = K,
      r = r_id, 
      s = s_id, 
      q = max_flow
    )

    C <- diag(K)
    C[r_idx, s_idx] <- 1 
    C <- C[, -r_idx, drop = FALSE] 
    
    F_matrix <- t(C) %*% F_matrix %*% C
    new_ids <- current_ids[-r_idx]
    rownames(F_matrix) <- colnames(F_matrix) <- new_ids

    F_matrix_history[[length(F_matrix_history) + 1]] <- list(
      F_matrix
    )

  }
  
  # --- 3. Generate Final Result Sets ---
  
  # --- 3a. Create unit_set ---
  unit_set <- data.frame(
    sourceunit = source_units,
    clusterid = as.character(source_units),
    destinationunit = NA,
    g = NA,
    round = NA,
    stringsAsFactors = FALSE
  )
  
  for (merge in rev(merge_history)) {
    is_merged_unit <- unit_set$clusterid == merge$r
    unit_set$clusterid[is_merged_unit] <- merge$s
    
    is_source_unit <- unit_set$sourceunit == merge$r
    unit_set$destinationunit[is_source_unit & is.na(unit_set$destinationunit)] <- merge$s
    unit_set$g[is_source_unit & is.na(unit_set$g)] <- merge$q
    unit_set$round[is_source_unit & is.na(unit_set$round)] <- merge$round
  }
  
  unit_set$core <- ifelse(is.na(unit_set$g), 1, 0)
  
  # --- 3b. Create cluster_set ---
  final_clusters <- rownames(F_matrix)
  
  row_flows <- rowSums(F_matrix)
  internal <- diag(F_matrix)
  
  internal_relative <- ifelse(row_flows == 0, 0, internal / row_flows)
  
  F_matrix_history <- lapply(F_matrix_history, \(x) x[[1]])
  names(F_matrix_history) <- max(unit_set$round,na.rm=TRUE):min(unit_set$round,na.rm=TRUE)

  cluster_set <- data.frame(
    clusterid = final_clusters,
    internal = internal,
    rowflows = row_flows,
    internal_relative = internal_relative,
    stringsAsFactors = FALSE
  )
  
  if(nrow(cluster_set) > 0) {
      cluster_set$La = mean(internal_relative)
      cluster_set$Lw = weighted.mean(internal_relative, w = row_flows)
      cluster_set$Lm = min(internal_relative)
      cluster_set$N = sum(row_flows)
  }

  if(save_k == TRUE) {
    return(list(unit_set = unit_set, cluster_set = cluster_set,
            F_matrix=F_matrix, F_matrix_history = F_matrix_history))
  } else {
    return(list(unit_set = unit_set, cluster_set = cluster_set,
            F_matrix=F_matrix))
  }
}

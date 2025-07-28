

flowbca_modularity <- function(unit_set, F_matrix_history){
 F_matrix_1 <- F_matrix_history[[1]]
 cluster_set <- cluster_set(unit_set)
 M_matrix <- CalcModMatrix(F_matrix_1)
 v_modularity <- unlist(lapply(cluster_set, \(x) calculate_modularity(F_matrix_1, M_matrix, x)))
 modularity <- data.frame(round=names(v_modularity), modularity=v_modularity)
 return(modularity)
}

cluster_set <- function(unit_set){
  unit_set <- as.data.frame(bca_result$unit_set)
  round_set <- unit_set[!is.na(unit_set$round), ]
  round_set <- round_set[order(round_set$round, decreasing = TRUE), ]

  cluster_set <- list()
  cluster_set[[1]] <- data.frame(sourceunit=unit_set$sourceunit,
                      clusterid=unit_set$sourceunit)
  rd <- round_set$round
  for(i in 1:length(rd)) {
    temp <- round_set[round_set$round==rd[i],]
    cluster_set[[i+1]] <- cluster_set[[i]]  
    cluster_set[[i+1]]$clusterid[cluster_set[[i+1]]$clusterid == temp$sourceunit] <- temp$destinationunit
  }
  names(cluster_set) <- c(max(rd+1),rd)
  return(cluster_set)
}

CalcModMatrix <- function(F_matrix_1, is_directed = FALSE, modularity_resolution = 1.0){
  mat <- F_matrix_1
  out_degree <- rowSums(mat)
  in_degree <- colSums(mat)
  total_weight <- sum(mat)
  E_degree <- outer(out_degree,in_degree) / total_weight
  m_modularity_matrix <- mat - modularity_resolution * E_degree
  return(m_modularity_matrix)
}

calculate_modularity <- function(F_matrix_1, M_matrix, communities){
  same_community_matrix <- outer(communities$clusterid, communities$clusterid, "==")
  modularity <- sum(M_matrix[same_community_matrix]) / sum(F_matrix_1)
  return(modularity)
}



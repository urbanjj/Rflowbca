# # (예시) 병합 로그 데이터프레임 df 준비
# df <- unit_set[order(unit_set$round), ]  # round 기준 정렬

# # 모든 노드 이름(leaf) 추출
# all_units <- unique(c(df$sourceunit))
# n_leaf <- length(all_units)
# labels <- all_units
# names(labels) <- 1:n_leaf

# # unit 이름 → 인덱스 맵핑
# unit2idx <- setNames(seq_len(n_leaf), all_units)
# current_idx <- n_leaf + 1

# # 각 클러스터가 현재 어떤 인덱스를 가지는지(업데이트)
# cluster_map <- unit2idx

# # 결과 merge/height 저장
# merge_mat <- matrix(NA, nrow=nrow(df), ncol=2)
# height_vec <- numeric(nrow(df))

# for(i in seq_len(nrow(df))) {
#   src <- as.character(df$sourceunit[i])
#   dst <- as.character(df$destinationunit[i])
  
#   # merge 행렬에 인덱스 기록 (leaf면 음수)
#   idx1 <- cluster_map[src]
#   idx2 <- cluster_map[dst]
#   merge_mat[i,] <- c(ifelse(idx1 <= n_leaf, -idx1, idx1 - n_leaf), 
#                      ifelse(idx2 <= n_leaf, -idx2, idx2 - n_leaf))
  
#   # height: round 또는 g 값 등
#   height_vec[i] <- df$round[i]
  
#   # 새로 만들어진 클러스터는 destinationunit의 이름으로 업데이트
#   cluster_map[dst] <- current_idx
#   # 사용된 sourceunit은 제거 (이미 합쳐졌으므로)
#   cluster_map[src] <- NA
#   current_idx <- current_idx + 1
# }

# # hclust 객체 만들기
# hc <- list(
#   merge = merge_mat,
#   height = height_vec,
#   order = 1:n_leaf,       # 실제 leaf 순서로 바꾸려면 order.dendrogram 등 활용
#   labels = labels,
#   method = "custom",
#   call = match.call()
# )
# class(hc) <- "hclust"

# # 그리기
# plot(hc, main="Custom hclust from merge log")

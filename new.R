## data
c1 <- read.csv('./c1.csv')[,-1]
c1$O_AD_CD <- paste0('X',c1$O_AD_CD)

# flowbca, build_hierarchy, build_cluster_tree, list_to_dendrogram, get_parent_node 함수는
# flowbcaR 패키지 내에 정의된 함수로 가정합니다.
c2 <- flowbca(c1)

d0 <- build_hierarchy(bca_result$unit_set)
d1 <- build_cluster_tree(d0)
aa <- hierarchy_cluster(d1)
flowbca_map(bca_result$unit_set, KR_SiGun, join_col=c('sourceunit'='SiGun_NM'))

c1 <- flowbca_map(aa[[1]]$hierarchy_1, KR_SiGun,
            join_col=c('sourceunit'='SiGun_NM'),filenm='tesf1.png')
c2 <- flowbca_map(aa[[1]]$hierarchy_2, KR_SiGun,
              join_col=c('sourceunit'='SiGun_NM'),filenm='tesf2.png')

unit_set <- aa[[1]]$hierarchy_2

d2 <- list_to_dendrogram(d1)

hierarchy_cluster <- function(cluster_tree){
  # list_to_dendrogram 함수는 패키지 환경에 정의된 것으로 가정합니다.

  flowbca_dendrogram <- list_to_dendrogram(cluster_tree)
  h <- attributes(flowbca_dendrogram)$height
  h_nm <- paste0('hierarchy_',h:1)

  h_list <- list()
  for(i in 1:(h-1)) {
    temp <- cut(flowbca_dendrogram, h=i)
    lower_list <- lapply(temp$lower, labels)
    names(lower_list) <- lapply(temp$lower, function(x) attr(x, "label"))
    df <- data.frame(
          sourceunit = unlist(lower_list),
          h_cl = rep(names(lower_list), times = sapply(lower_list, length)),
          row.names = NULL,
          stringsAsFactors = FALSE)
    h_list[[i]] <- df
  }
  names(h_list) <- h_nm[-1]

  cluster_info <- list()
  core_info <- list()
  for(i in 1:(h-1)){
    temp_df <- h_list[[i]]
    
    # 멤버가 1개 이상인 클러스터(코어)를 찾습니다.
    h_cl_counts <- table(temp_df$h_cl)
    multi_member_clusters <- names(h_cl_counts[h_cl_counts > 1])
    df <- temp_df[temp_df$h_cl %in% multi_member_clusters, ]

    # get_parent_node 함수는 패키지 환경에 정의된 것으로 가정합니다.
    # 전역 변수 d1 대신 인자로 받은 cluster_tree를 사용하도록 수정했습니다.
    p1 <- unlist(lapply(df$h_cl, function(x) get_parent_node(cluster_tree, x)))
    cluster_info[[i]] <- cbind(df, upper_h_cl=NA)
    cluster_info[[i]]$upper_h_cl <- p1

    # 코어 클러스터에 대한 조회 테이블을 생성합니다.
    core_counts_df <- as.data.frame(table(temp_df$h_cl), stringsAsFactors = FALSE)
    names(core_counts_df) <- c("sourceunit", "core_N")
    df_core <- core_counts_df[core_counts_df$core_N > 1, ]
    if(nrow(df_core) > 0) {
      df_core$core <- 1
    }
    core_info[[i]] <- df_core
  }
  names(cluster_info) <- h_nm[-1]
  names(core_info) <- h_nm[-1]
  
  # 원본 코드의 동작을 그대로 재현합니다.
  # h_list의 'sourceunit'(단위 ID)와 core_info의 'sourceunit'(클러스터 이름)을 기준으로 병합합니다.
  # 의도와 다를 수 있으나, 원본 코드의 동작을 유지하기 위함입니다.
  unit_set <- Map(function(x,y) {
                                 rl <- merge(x, y, by = 'sourceunit', all.x = TRUE)
                                 # 병합 과정에서 생긴 NA를 0으로 바꿉니다.
                                 rl[is.na(rl)] <- 0
                                 return(rl)
                                 }, h_list, core_info)
  
  hierarchy_cluster <- list(unit_set, cluster_info)
  names(hierarchy_cluster) <- c('unit_set','cluster_info')
  return(hierarchy_cluster)
}

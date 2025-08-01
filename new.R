## data
c1 <- read.csv('./c1.csv')[,-1]
c1$O_AD_CD <- paste0('X',c1$O_AD_CD)
c2 <- flowbca(c1)

d0 <- build_hierarchy(c2$unit_set)

d0 <- build_hierarchy(bca_result$unit_set)
d1 <- build_cluster_tree(d0)
d2 <- list_to_dendrogram(d1)

hierarchy_cluster <- function(flowbca_dendrogram){}

flowbca_dendrogram <- d2
h <- attributes(flowbca_dendrogram)$height
h_nm <- paste0('hierarchy_',h:1)

h_list <- list()
for(i in 1:(h-1)) {
  temp <- cut(flowbca_dendrogram, h=i)
  lower_list <- lapply(temp$lower, labels)
  names(lower_list) <- lapply(temp$lower, \(x) attr(x, "label")) 
  df <- data.frame(
        sourceunit = unlist(lower_list),
        h_cl = rep(names(lower_list), times = lapply(lower_list, length)),
        row.names = NULL)
  h_list[[i]] <- df
}
names(h_list) <- h_nm[-1]

m_list <- list()
n_list <- list()
for(i in 1:(h-1)){
  temp <- as.data.table(h_list[[i]])
  df <- temp[h_cl %in% temp[,.N,by=.(h_cl)][N!=1]$h_cl]
  p1 <- unlist(lapply(df$h_cl, \(x) get_parent_node(d1, x)))
  m_list[[i]] <- cbind(df, upper_h_cl=p1)
  n_list[[i]] <- temp[,.N,by=.(h_cl)][N!=1]
}
names(m_list) <- h_nm[-1]
names(n_list) <- h_nm[-1]


m0 <- m_list[[2]]
n0 <- n_list[[2]]

KR_SiGun
r1 <- merge(KR_SiGun, m0, by.x='SiGun_NM', by.y='sourceunit', all.y=TRUE)

library(dplyr)
library(tidyverse)
r2 <- r1 %>% group_by(h_cl) %>% summarise
plot(KR_SiGun$geom)
plot(r2, add=TRUE)























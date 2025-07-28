library(tidyverse)
library(fastverse)
library(igraph)

## igraph
g <- make_graph("Zachary")

Z1 <- igraph::as_data_frame(Z0, what = "edges") %>%
  mutate(from=paste0('z',sprintf("%02d",from)),to=paste0('z',sprintf("%02d",to)),weight=1)

Z2 <- expand_grid(from=paste0('z',sprintf("%02d",1:34)),
                  to=paste0('z',sprintf("%02d",1:34))) %>%
  left_join(Z1, by=c('from','to')) %>%
  mutate(weight=ifelse(is.na(weight),0,weight))



F_matrix <- dcast(as.data.table(Z2), from~to, value.var='weight', fill=0)
F_matrix_1 <- F_matrix[,-1]
rownames(F_matrix_1) <- colnames(F_matrix_1)

Zachary_membership <- igraph::cluster_infomap(g)
membership(Zachary_membership)
Zachary_membership$modularity

qq <- modularity_matrix(g)
(M_matrix - qq) %>% write_clip


u0 <- data.frame(sourceunit=paste0('z',sprintf("%02d",1:34)),
        clusterid=membership(Zachary_membership))

M_matrix <- CalcModMatrix(adj_matrix, is_directed = TRUE)
zz <- calculate_modularity(F_matrix_1, M_matrix, u0, is_directed = TRUE)
zz


### 
library(clipr)
library(igraph)
g <- make_graph("Zachary")
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE, )
write_clip(adj_matrix)

g_dir <- as_directed(g, mode = "mutual")  # 또는 "arbitrary"
is_directed(g_dir)  

el <- as_edgelist(g, names = FALSE)
n <- vcount(g)
m <- matrix(0, nrow = n, ncol = n)
for (i in 1:nrow(el)) {
  m[el[i,1], el[i,2]] <- 1
}
isSymmetric(m)

is_directed(g)

i_M_matrix <- modularity_matrix(g)
M_matrix <- CalcModMatrix(adj_matrix, is_directed = TRUE)
sum(i_M_matrix - M_matrix)

mem <- cluster_infomap(g)

mem$modularity
calculate_modularity(adj_matrix, M_matrix, mem$membership)
modularity(g, mem$membership, directed=FALSE)

mem2 <- cluster_louvain(g)
mem2$modularity
calculate_modularity(adj_matrix, M_matrix, mem2$membership)

mem3 <- cluster_edge_betweenness(g)
mem3$modularity
calculate_modularity(adj_matrix, M_matrix, mem3$membership)



M_matrix <- CalcModMatrix(adj_matrix, is_directed = FALSE)
calculate_modularity(adj_matrix, M_matrix, mem$membership, is_directed = FALSE)
modularity(g, mem$membership, directed=FALSE)


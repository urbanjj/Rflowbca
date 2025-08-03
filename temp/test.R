library(tidyverse)
library(fastverse)
library(igraph)
library(clipr)


g <- make_graph("Zachary")
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE, )
write_clip(adj_matrix)

i_M_matrix <- modularity_matrix(g, resolution = 4)
M_matrix <- CalcModMatrix(adj_matrix, is_directed = TRUE, resolution = 4)
sum(i_M_matrix - M_matrix)

mem <- cluster_infomap(g)
mem$modularity
calculate_modularity(adj_matrix, mem$membership)
calculate_modularity(adj_matrix, mem$membership,4)
modularity(g, mem$membership, directed=FALSE,resolution=1)
modularity(g, mem$membership, directed=FALSE,resolution=4)


mem2 <- cluster_louvain(g)
mem2$modularity
calculate_modularity(adj_matrix, mem2$membership)

mem3 <- cluster_edge_betweenness(g)
mem3$modularity
calculate_modularity(adj_matrix, mem3$membership)

adj_matrix
M_matrix <- CalcModMatrix(adj_matrix, is_directed = FALSE)
calculate_modularity(adj_matrix, mem$membership, is_directed = FALSE)
modularity(g, mem$membership, directed=FALSE)

adj_matrix2 <- adj_matrix
rownames(adj_matrix2) <- paste0('X',c(1:34))
colnames(adj_matrix2) <- paste0('X',c(1:34))

CalcModMatrix(adj_matrix, is_directed = FALSE)
mem_v <- mem$membership
names(mem_v) <- paste0('X',c(1:34))
calculate_modularity(adj_matrix2, mem_v, is_directed = FALSE)

summary(rownames(adj_matrix2) == colnames(adj_matrix2)) == TRUE

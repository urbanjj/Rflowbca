library(tidyverse)
library(fastverse)
library(igraph)
library(clipr)


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


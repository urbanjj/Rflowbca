nested_list <- list(Nation=c(Daegu=list(AA0$Daegu),list(Busan=AA0$Busan)))

nested_list2 <- list(Nation=AA0)

r1 <- list_to_dendrogram(nested_list2)

labels(r1)

str(r1)

r2 <- cut(r1, h=2)
r2$upper %>% str()
sapply(r2$upper, labels) 

## name 확인
sapply(r2$lower, \(x) attr(x, "label")) 

sapply(r2$lower, \(x) attr(x, "label")) 
sapply(r2$lower, labels) 



r2$lower %>% str()

plot(r1)

## data
c1 <- read.csv('./c1.csv')[,-1]
c1$O_AD_CD <- paste0('X',c1$O_AD_CD)
c2 <- flowbca(c1)

d0 <- build_hierarchy(c2$unit_set)
d1 <- build_cluster_tree(d0)
d2 <- list_to_dendrogram(list(nation=d1))

hierarchy_cluster <- function(flowbca_dendrogram){}
flowbca_dendrogram <- d2
h <- attributes(flowbca_dendrogram)$height

h_list <- list()
for(i in 1:(h-1)) {
  temp <- cut(flowbca_dendrogram, h=i)
  lower_list <- sapply(temp$lower, labels)
  names(lower_list) <- sapply(temp$lower, \(x) attr(x, "label")) 
  df <- data.frame(
        sourceunit = unlist(lower_list),
        h_cl = rep(names(lower_list), times = sapply(lower_list, length)),
        row.names = NULL)
  h_list[[i]] <- df
}

plot(d2)

m0 <- as.data.table(h_list[[4]])
m0[,.N,by=.(h_cl)][N!=1]
m1 <- m0[h_cl %in% m0[,.N,by=.(h_cl)][N!=1]$h_cl]
p1 <- unlist(lapply(m1$h_cl, \(x) get_parent_node(d1, x)))
m1[,upper_h_cl:=p1]


get_parent_node(d1, "Gumi")


proune()
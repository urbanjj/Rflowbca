#' Find Hierarchy Path using Base R (Internal Helper)
#'
#' Traverses a parent-child map to find the full hierarchy path for a single node.
#' This is an internal helper function implemented using only base R.
#'
#' @param start_node The node to start the traversal from.
#' @param cluster_map A named vector where names are children and values are parents.
#' @return A string representing the full path.
find_hierarchy <- function(start_node, cluster_map) {
  path <- c(start_node)
  current_node <- start_node
  visited <- c(start_node)

  while (current_node %in% names(cluster_map)) {
    parent_node <- cluster_map[[current_node]]

    # Stop if the parent is NA or empty
    if (is.na(parent_node)) {
      break
    }

    # Check for circular references
    if (parent_node %in% visited) {
      warning(paste("Circular reference detected:", paste(c(visited, parent_node), collapse=" -> ")))
      path <- c(path, "...") # Indicate a circular path
      break
    }

    path <- c(path, parent_node)
    visited <- c(visited, parent_node)
    current_node <- parent_node
  }

  return(paste(rev(path), collapse = "/"))
}

#' Build and Add Hierarchy Path Column using Base R
#'
#' This function constructs hierarchy paths from a data frame of parent-child
#' relationships using only base R functions. It is a replacement for the
#' original `find_hierarchy_db`.
#'
#' @param data A data frame, typically the `unit_set` from a `flowbca` result.
#' @param child_col A string, the name of the column with child units (defaults to "sourceunit").
#' @param parent_col A string, the name of the column with parent units (defaults to "destinationunit").
#' @return The input data frame with a new `hierarchy` column.
#' @importFrom stats setNames
#' @export
build_hierarchy <- function(data, child_col = "sourceunit", parent_col = "destinationunit") {
  # Ensure input is a data frame
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  # Create the parent-child map from the data frame
  valid_links <- data[!is.na(data[[parent_col]]), ]
  cluster_map <- setNames(valid_links[[parent_col]], valid_links[[child_col]])
  clusterid <- unique(data$clusterid)

  # Use sapply to apply the find_hierarchy function to each child unit
  data$hierarchy <- sapply(data[[child_col]], find_hierarchy, cluster_map = cluster_map)
  data$h_level <- sapply(data$hierarchy, \(x) nchar(x)-nchar(gsub('/','',x)))+1

  # Calculate the number of children for each node
  child_counts <- table(data[[parent_col]])
  data$child_count <- as.integer(child_counts[as.character(data[[child_col]])])
  data$child_count[is.na(data$child_count)] <- 0
  parentid <- data$sourceunit[data$child_count > 0]
  data$h_parent <- sapply(data$sourceunit, \(x) ifelse(x %in% parentid,1,0))

  return(data)
}

#' Build a Nested List Representing the Cluster Hierarchy
#'
#' This function converts the flat parent-child relationship in `unit_set` into a
#' nested list (tree) structure. This is useful for hierarchical analysis or
#' visualization with tree-based packages.
#'
#' @param unit_set A data frame from a `flowbca` result, containing at least
#'   `sourceunit` and `destinationunit` columns defining the parent-child links.
#' @return A named, nested list representing the cluster hierarchy. Each name is a
#'   node, and its value is a list of its children. Leaf nodes are empty lists.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'bca_result' is the output from flowbca()
#'   tree <- build_cluster_tree(bca_result$unit_set)
#'   # You can inspect the tree structure
#'   str(tree)
#'   # Access a specific sub-tree
#'   seoul_cluster <- tree$Seoul
#' }

build_cluster_tree <- function(unit_set) {

  # 1. Create a simple, reliable adjacency list (parent -> vector of children)
  valid_links <- unit_set[!is.na(unit_set$destinationunit), ]
  if (nrow(valid_links) == 0) {
    # Handle case with no merges by returning a list of all units as empty lists
    all_units <- unit_set$sourceunit[!is.na(unit_set$sourceunit)]
    return(setNames(lapply(all_units, function(x) list()), all_units))
  }
  adjacency_list <- split(as.character(valid_links$sourceunit), as.character(valid_links$destinationunit))

  # 2. Define the recursive function to build the nested list from the adjacency list
  build_recursive <- function(node_name) {
    children <- adjacency_list[[node_name]]
    if (is.null(children) || length(children) == 0) {
      return(list())  # leaf node: do not include itself as child
    }
    children <- sort(children)
    child_list <- setNames(lapply(children, build_recursive), children)
    # If the node has at least one child, add itself as a child
    child_list[[node_name]] <- list()
    # Sort by name (optional)
    child_list[sort(names(child_list))]
  }

  # 3. Identify the top-level nodes (roots of the trees)
  # These are nodes that are not children of any other node.
  all_units <- unique(c(unit_set$sourceunit, unit_set$destinationunit))
  all_units <- all_units[!is.na(all_units)]
  merged_units <- unit_set$sourceunit[!is.na(unit_set$destinationunit)]
  top_level_nodes <- sort(setdiff(all_units, merged_units))

  # 4. Build the final tree by calling the recursive function on each top-level node
  final_tree <- setNames(lapply(top_level_nodes, build_recursive), top_level_nodes)
  return(final_tree)
}

#' Extract and Organize Clusters at a Specific Hierarchy Level
#'
#' This function identifies clusters at a given hierarchy level and organizes them
#' based on their parent-child relationships. For a specified level, it returns
#' a list where each element corresponds to a cluster.
#'
#' If a node at the specified level is a child of a node from the level above,
#' the parent node's cluster is returned with the child node's members excluded.
#' The child node itself is returned as a separate cluster including all its members.
#'
#' @param unit_set A data frame from a `flowbca` result, which must contain
#'   `sourceunit` and `destinationunit` columns to define the hierarchy.
#' @param hierarchy An integer specifying the hierarchy level to analyze (must be >= 1).
#'   Defaults to 1.
#' @return A named list of character vectors. Each name is the root of a cluster,
#'   and the value is a vector of units belonging to that cluster, adjusted for
#'   sub-cluster relationships.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'bca_result' is the output from flowbca()
#'   # Get clusters at level 2
#'   clusters_level_2 <- cluster_hierarchy(bca_result$unit_set, hierarchy = 2L)
#'
#'   # This will return a list where, for example, a parent cluster like 'Seoul'
#'   # might appear, but with members of its sub-cluster at level 2 (e.g., 'Incheon')
#'   # removed. 'Incheon' would then be a separate item in the list with its own members.
#' }
cluster_hierarchy <- function(unit_set, hierarchy = 2L) {
  hierarchy <- as.integer(hierarchy)
  if (!is.data.frame(unit_set) || !all(c("h_parent", "h_level") %in% names(unit_set))) {
    stop("unit_set must be a data frame with 'sourceunit' and 'destinationunit' columns.")
  }
  if (!is.integer(hierarchy) || hierarchy < 1 || hierarchy > max(unit_set$h_level)-1) {
    stop("hierarchy must be an integer between 1 and max of h_level.")
  }

  # 1. Build hierarchy data and identify nodes at the target and parent levels
  cluster_tree <- build_cluster_tree(unit_set)
  
  if(hierarchy == 1){
    hchy_1_nm <- unit_set$sourceunit[unit_set$h_parent==1 & unit_set$h_level==hierarchy]
    hchy_1_des <- lapply(as.list(hchy_1_nm), \(x) get_descendants_self(x, cluster_tree))
    names(hchy_1_des) <- hchy_1_nm
    return(hchy_1_des)
  }
  
  if(hierarchy > 1){
    hchy_up_nm <- unit_set$sourceunit[unit_set$h_parent==1 & unit_set$h_level==(hierarchy-1)]
    hchy_down_nm <- unit_set$sourceunit[unit_set$h_parent==1 & unit_set$h_level==hierarchy]
    
    hchy_up_parent <- lapply(as.list(hchy_up_nm), \(x) get_parent_node(cluster_tree,x))
    hchy_up_des <- Map(\(x,y) get_descendants_self(x,cluster_tree[[y]]), as.list(hchy_up_nm), hchy_up_parent)
    hchy_up_des_exclude <- lapply(as.list(hchy_up_nm),\(x) c(x,exclude_descendants(cluster_tree[[x]], hchy_down_nm)))
    names(hchy_up_des_exclude) <- hchy_up_nm

    hchy_down_parent <- lapply(as.list(hchy_down_nm), \(x) get_parent_node(cluster_tree,x))
    hchy_down_des <- Map(\(x,y) get_descendants_self(x,cluster_tree[[y]]), as.list(hchy_down_nm), hchy_down_parent)
    names(hchy_down_des) <- hchy_down_nm
    
    cluster_hierarchy <- list('root_exclude'=hchy_up_des_exclude, 'sub'=hchy_down_des)
    return(cluster_hierarchy)
  }
}


#' Get All Descendants from a Nested List (Helper)
#'
#' Recursively traverses a nested list (tree) structure to find all child nodes
#' (descendants) starting from a given node.
#'
#' @param node A list representing the current node in the hierarchy tree.
#' @return A character vector containing the names of all descendant nodes.
#' @noRd
get_all_descendants_nested <- function(node) {
  if (length(node) == 0) {
    return(character(0))
  }
  result <- character(0)
  for (child_name in names(node)) {
    result <- c(result, child_name, get_all_descendants_nested(node[[child_name]]))
  }
  result
}


#' Get a Node and All Its Descendants (Helper)
#'
#' A wrapper function that returns the starting node itself along with all of its
#' descendants by calling `get_all_descendants_nested`.
#'
#' @param node_name The name of the starting node (character).
#' @param node The full hierarchy tree (nested list) to search within.
#' @return A character vector containing the starting node and all its descendants.
#' @noRd
get_descendants_self <- function(node_name, node) {
  descendants <- get_all_descendants_nested(node[[node_name]])
  c(node_name,descendants)
}

#' Exclude Specific Descendant Trees from a Parent Node (Helper)
#'
#' From a given parent node's list of all descendants, this function removes
#' one or more specified child nodes and all of their own descendants.
#'
#' @param parent_node The parent node's sub-tree (a list).
#' @param exclude_nodes A character vector of child nodes to exclude.
#' @return A character vector of descendants of the parent, with the specified
#'   sub-trees removed.
#' @noRd
exclude_descendants <- function(parent_node, exclude_nodes) {
  all_descendants <- get_all_descendants_nested(parent_node)
  exclude_all <- unlist(lapply(exclude_nodes, function(n) get_descendants_self(n, parent_node)))
  setdiff(all_descendants, exclude_all)
}


#' Find the Parent of a Given Node (Helper)
#'
#' Recursively searches the entire cluster tree to find the immediate parent
#' of a specified target node.
#'
#' @param node The current node or sub-tree (list) to search in.
#' @param target_name The name of the node whose parent is to be found.
#' @param parent The name of the parent in the current recursion level (internal use).
#' @return The name of the parent node (character) or NULL if not found.
#' @noRd
get_parent_node <- function(node, target_name, parent = NULL) {
  for (name in names(node)) {
    if (name == target_name) {
      return(parent)
    }
    child <- node[[name]]
    if (is.list(child) && length(child) > 0) {
      found <- get_parent_node(child, target_name, name)
      if (!is.null(found)) return(found)
    }
  }
  return(NULL)
}

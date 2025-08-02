## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(flowbcaR)
library(sf)

## ----workflow-function, echo=FALSE, include = FALSE, eval =FALSE--------------
#  library(DiagrammeR)
#  library(DiagrammeRsvg)
#  
#  diagram <- grViz("
#  digraph flowbcaR_flow {
#  
#    graph [rankdir = LR]
#  
#    subgraph cluster_input {
#      label = '1. Input Data';
#      style = dashed;
#      A [label = 'flow_input\\n(OD data.frame)', shape = folder]
#      B [label = 'KR_SiGun\\n(sf object)', shape = folder]
#    }
#  
#    subgraph cluster_core {
#      label = '2. Core Clustering';
#      C [label = 'flowbca', shape = box, style = filled, fillcolor = lightblue]
#    }
#  
#    subgraph cluster_analysis {
#      label = '3. Analysis & Visualization';
#      D [label = 'build_hierarchy', shape = box]
#      E [label = 'flowbca_stat', shape = box]
#      F [label = 'flowbca_plot', shape = box]
#      G [label = 'flowbca_modularity', shape = box]
#      H [label = 'flowbca_gis', shape = box]
#      I [label = 'flowbca_ani', shape = box]
#    }
#  
#    subgraph cluster_diagnosis {
#      label = '4. All-in-One Diagnosis';
#      J [label = 'flowbca_diagnosis', shape = box, style = filled, fillcolor = lightgray]
#    }
#  
#    # Outputs
#    C_out [label = 'bca_result\\n(list)', shape = note, color = blue]
#    E_out [label = 'stat_data', shape = note, color = blue]
#    H_out [label = 'gis_layers', shape = note, color = blue]
#  
#    # Edges
#    A -> C
#    C -> C_out
#  
#    C_out -> D [label = 'unit_set']
#    C_out -> E [label = 'F_matrix_history']
#    E -> E_out
#    E_out -> F
#  
#    C_out -> G [label = 'unit_set &\\nF_matrix_history']
#  
#    C_out -> H [label = 'unit_set']
#    B -> H
#    H -> H_out
#    H_out -> I [label = 'gis_layers']
#    C_out -> I [label = 'unit_set']
#  
#    A -> J
#  }
#  ")
#  diagram_svg <- export_svg(diagram)
#  writeLines(diagram_svg, "workflow.svg")

## ----workflow_plot, echo=FALSE, out.width="100%", fig.align = 'center'--------
knitr::include_graphics("../man/figures/workflow.svg")

## ----prepare-data, echo=TRUE--------------------------------------------------
library(flowbcaR)
library(sf)

# Load the sample datasets
data(OD_SiGun)
data(KR_SiGun)

# Prepare the flow data for the algorithm
# The first column must be the source unit ID, and the destination column
# names must match the source unit IDs.
# We remove the first column ("SiGun_CD") to meet this requirement.
flow_input <- OD_SiGun[, -1]
colnames(flow_input) <- c('SiGun_NM',flow_input[,1])

# Check the prepared data
print("Prepared OD Data for flowbca:")
print(flow_input[1:5, 1:6])
print(head(KR_SiGun))

## ----flowbca-example, cache=TRUE----------------------------------------------
# Run clustering until the minimum internal relative flow (lm) is at least 0.1 (10%).
# save_k=TRUE is important to retain data for subsequent analysis functions.
bca_result <- flowbca(flow_input, lm = 0.1, save_k = TRUE)

# The result is a list containing:
# 1. unit_set: Details of cluster assignment for each unit.
# 2. cluster_set: Statistics for the final clusters.
# 3. F_matrix: The final aggregated flow matrix.
# 4. F_matrix_history: A list of matrices from each clustering round.
# 5. C_matrix_history: A list of transformation matrices.
str(bca_result$unit_set, 3)
str(bca_result$cluster_set, 3)

## ----hierarchy-example--------------------------------------------------------
# Build the hierarchy path from the result
hierarchy_data <- build_hierarchy(bca_result$unit_set)

# View the hierarchy for a few units
head(hierarchy_data[, c("sourceunit", "clusterid", "hierarchy", "h_level")])

## ----data-tree-example, message=FALSE-----------------------------------------
library(data.tree)
# Create a pathString for data.tree
hierarchy_data$pathString <- paste("Korea", hierarchy_data$hierarchy, sep="/")
tree <- as.Node(hierarchy_data, pathName="pathString")
# Print the top levels of the tree
print(tree, "level", limit=10)
# Explore the sub-tree structure for Busan
print(tree$Busan)

## ----stat-plot-example, out.width="90%", fig.width=10, fig.height=5, fig.align = 'center'----
# Calculate internal flow statistics from the matrix history
stat_data <- flowbca_stat(bca_result$F_matrix_history)
str(stat_data)
# Plot the statistics
# The x-axis represents the round (number of clusters + 1)
# The plot is interactive when the number of points is small (<= 20)
flowbca_plot(stat_data)

## ----modularity-example, out.width="60%", fig.width=8,fig.height=6, fig.align = 'center'----
# Calculate modularity for each round
modularity_data <- flowbca_modularity(bca_result$unit_set, bca_result$F_matrix_history)

# Plot modularity over rounds
plot(modularity_data$round, modularity_data$modularity, type='l',
        xlab="Round", ylab="Modularity", main="Modularity vs. Clustering Step",
        xlim=rev(range(as.integer(modularity_data$round))))

## ----diagnosis-example, out.width="90%", fig.width=8,fig.height=7, fig.align = 'center'----
# Perform diagnosis for a directed graph
# This function calls flowbca internally, so it may take a moment
diagnosis_stat <- flowbca_diagnosis(flow_input, is_directed = TRUE)
# The function returns the statistics used for plotting
str(diagnosis_stat)


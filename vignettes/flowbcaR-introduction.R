## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(flowbcaR)
library(sf)

## ----workflow-diagram-png, echo=FALSE, out.width="100%"-----------------------
knitr::include_graphics("../man/figures/workflow-diagram.png")

## ----prepare-data, echo=TRUE--------------------------------------------------
# Load the sample datasets
data(OD_SiGun)
data(KR_SiGun)

# Prepare the flow data: use SiGun_NM as the ID
# Remove the numeric code column and set column names correctly.
flow_input <- OD_SiGun[, -1]
# The first column name must match the pattern of the destination columns.
# Here we rename it to be descriptive, but it's not used by the function.
colnames(flow_input) <- c('SiGun_NM', OD_SiGun$SiGun_NM)


# Check the prepared data
print("Prepared OD Data for flowbca:")
print(flow_input[1:5, 1:6])


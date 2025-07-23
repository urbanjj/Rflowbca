#' Origin-Destination Flow Data for Si-Gun in South Korea
#'
#' A dataset containing origin-destination (OD) flow data between Si (cities)
#' and Gun (counties) in the Republic of Korea. It is structured as a data frame
#' where each row represents an origin, and columns represent destinations.
#'
#' @name OD_SiGun
#' @docType data
#' @format A \code{data.frame} with 159 observations and 161 variables:
#' \describe{
#' \item{SiGun_CD}{A character string representing the unique administrative code of the origin unit.}
#' \item{SiGun_NM}{A character string in UTF-8 representing the Korean name of the origin unit.}
#' \item{...}{The subsequent 159 columns are named with destination \code{sourceunit} codes (e.g., "11", "21"). The numeric values in these columns represent the flow from the origin (row) to the destination (column).}
#' }
#'
#' @source Fundamental transport demand data provided by the Korea Transport Institute. It contains total passenger origin-destination (O/D) trip volumes across all trip purposes.
#'
#' @keywords datasets flow od korea
#'
#' @examples
#' # To load the data from the package:
#' data(OD_SiGun)
#' # View the first few rows and columns
#' OD_SiGun[1:5, 1:5]
"OD_SiGun"

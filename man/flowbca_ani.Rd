% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/animation.R
\name{flowbca_ani}
\alias{flowbca_ani}
\title{Create a GIF Animation of the Clustering Process}
\usage{
flowbca_ani(
  flowbca_gis_layer,
  unit_set,
  width = 1000,
  filenm = "flowbca_animation.gif",
  keep_frames = FALSE
)
}
\arguments{
\item{unit_set}{A data frame from the \code{flowbca} result, containing details
about the merges, including \code{round}, \code{sourceunit}, and \code{clusterid}.}

\item{width}{The width of the output GIF in pixels. The height is calculated
automatically to maintain the aspect ratio. Defaults to 1000.}

\item{filenm}{The file name for the output GIF. Defaults to \code{'flowbca_animation.gif'}.}

\item{keep_frames}{A logical value. If \code{TRUE}, the individual PNG frames used to
create the GIF are kept in a sub-directory. Defaults to \code{FALSE}.}

\item{flowbca_gis}{A named list of \code{sf} objects, typically the output from the
\code{flowbca_gis_layer} function.}
}
\value{
This function does not return a value but writes a GIF file to the specified
location.
}
\description{
Generates a GIF animation visualizing the step-by-step merging of spatial units
from a \code{flowbca} analysis. It creates a frame for each round of the clustering
process and combines them into a single animated GIF.
}
\note{
Requires the \code{sf} and \code{magick} packages. For faster and
higher-quality PNG generation, it is recommended to also install the \code{ragg} package.
}

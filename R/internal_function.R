#' Remove Holes from a Single sf Geometry
#'
#' This function takes a single sf geometry (POLYGON or MULTIPOLYGON) and
#' removes any interior holes by reconstructing it using only the exterior ring(s).
#'
#' @param geom An `sf` geometry object (e.g., from `st_geometry()`).
#' @return An `sf` geometry object of the same type but without interior holes.
#' @noRd
remove_holes <- function(geom) {
  if (inherits(geom, "POLYGON")) {
    sf::st_polygon(list(geom[[1]]))
  } else if (inherits(geom, "MULTIPOLYGON")) {
    sf::st_multipolygon(lapply(geom, function(p) list(p[[1]])))
  } else {
    geom
  }
}

#' Remove Interior Holes from all Geometries in an sf Object
#'
#' This function iterates over all geometries within an `sf` object and removes
#' any interior holes (islands) from POLYGON and MULTIPOLYGON features. It is a
#' wrapper around the internal `remove_holes` function.
#'
#' @param sf_obj An `sf` object containing POLYGON or MULTIPOLYGON geometries.
#' @return A new `sf` object with holes removed from its geometries.
#' @export
remove_holes_sf <- function(sf_obj) {
  geom <- sf::st_geometry(sf_obj)
  new_geom <- sf::st_sfc(lapply(geom, remove_holes), crs = sf::st_crs(sf_obj))
  sf::st_set_geometry(sf_obj, new_geom)
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

#' @title Deprecated label spacing functions
#'
#' @description
#' `spaceLabels_new_new()` was deprecated in favor of `getLabelCoordinates()`.
#' `spaceLabels_newww()` and `spaceLabels_new()` were deprecated in favor of `spaceLabels()`.
#'
#' @name MyQuakyTools-deprecated
#' @keywords internal
NULL

#' @rdname MyQuakyTools-deprecated
#' @export
spaceLabels_new_new <- function(axis = "y-axis", max.distance = 6,
                                arrow.coordinates, label.coordinates, minDist) {

  # The warning: "version_number", "old_function()", "new_function()"
  lifecycle::deprecate_warn("1.1.0", "MyQuakyTools::spaceLabels_new_new()", "getLabelCoordinates()")

  # Pass arguments to the new function
  getLabelCoordinates(axis = axis, max.distance = max.distance,
                      arrow.coordinates = arrow.coordinates,
                      label.coordinates = label.coordinates, minDist = minDist)
}

#' @rdname MyQuakyTools-deprecated
#' @export
spaceLabels_newww <- function(labelCoords, minDist = 0.4) {

  lifecycle::deprecate_warn("1.1.0", "MyQuakyTools::spaceLabels_newww()", "spaceLabels()")

  spaceLabels(labelCoords = labelCoords, minDist = minDist)
}

#' @rdname MyQuakyTools-deprecated
#' @export
plotVolcalno <- function(...) {

  lifecycle::deprecate_warn("1.1.0", "MyQuakyTools::plotVolcalno()", "plotVolcano()")

  plotVolcano(...)
}

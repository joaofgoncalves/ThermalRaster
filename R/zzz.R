#' Package Attachment Welcome Message
#'
#' Displays a welcome message with the current package version when the package is loaded.
#'
#' @param libname The library path of the attached package (automatically provided by R).
#' @param pkgname The name of the package being attached (automatically provided by R).
#'
#' @keywords internal
#' @importFrom utils packageVersion
#' @noRd
#'
.onAttach <- function(libname, pkgname) {
  # Retrieve the package version
  version <- as.character(utils::packageVersion("ThermalRaster"))

  # Create a custom message with the package version
  message <- paste0("Welcome to ThermalRaster version: ", version, "!")

  # Display the startup message
  packageStartupMessage(message)
}

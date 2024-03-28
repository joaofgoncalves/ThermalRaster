
.onAttach <- function(libname, pkgname) {
  # Retrieve the package version
  version <- as.character(packageVersion("ThermalRaster"))

  # Create a custom message with the package version
  message <- paste0("Welcome to ThermalRaster version: ", version, "!")

  # Display the startup message
  packageStartupMessage(message)
}

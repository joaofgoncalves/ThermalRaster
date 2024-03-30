
#' Create a mock raster for testing
#'
#' Generates a mock raster object with specified dimensions and values for testing purposes.
#' This function is useful for creating sample raster data when developing or testing
#' spatial analysis workflows.
#'
#' @param values A vector of values to fill the raster. The length of this vector
#'               should be equal to `nrows` * `ncols`. Values are filled in row-wise.
#' @param nrows The number of rows in the raster.
#' @param ncols The number of columns in the raster.
#'
#' @return A `SpatRaster` object with the specified dimensions and values.
#'
#' @importFrom terra rast
#'
#' @export
#'
create_mock_raster <- function(values, nrows, ncols) {
  rast <- rast(nrows=nrows, ncols=ncols, vals=values)
  return(rast)
}

test_that("sobel_filter returns expected structure", {
  # For raster creation and operations# Create a simple 3x3 raster with a clear edge
  values <- c(255, 255, 255,
              0, 0, 0,
              255, 255, 255)
  mock_rast <- create_mock_raster(values, 3, 3)
  result <- sobel_filter(mock_rast)

  # Check if result is a terra Raster object
  expect_true(is(result, "SpatRaster"))

})


test_that("sobel_filter handles empty (NA) raster", {
  # Create an empty raster
  mock_rast <- create_mock_raster(rep(NA, 9), 3, 3)
  result <- sobel_filter(mock_rast)

  # Expect all values to be NA
  expect_true(all(is.na(values(result))))
})


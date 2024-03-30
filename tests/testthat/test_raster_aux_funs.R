

# Helper function to create a mock raster
create_mock_raster <- function(values) {
  rast <- rast(nrows=1, ncols=length(values), vals=values)
  return(rast)
}

test_that("remove_outliers correctly handles default thresholds", {
  values <- c(1, 2, 3, 4, 95, 96, 97, 98, 99, 100)
  mock_rast <- create_mock_raster(values)
  result <- remove_outliers(mock_rast)
  result_values <- values(result)[,1]

  # Expect values outside the 5th and 100th percentiles to be NA
  expected <- c(NA, 2, 3, 4, 95, 96, 97, 98, 99, 100)
  expect_equal(result_values, expected)
})
#
# test_that("remove_outliers works with custom thresholds", {
#   values <- c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#   mock_rast <- create_mock_raster(values)
#   result <- remove_outliers(mock_rast, pmin = 10, pmax = 90)
#   result_values <- values(result)
#
#   # Expect values outside the 10th and 90th percentiles to be NA
#   expected <- c(NA, 10, 20, 30, 40, 50, 60, 70, 80, 90, NA)
#   expect_equal(result_values, expected)
# })
#
# test_that("remove_outliers handles rasters with all NAs when thresholds are extreme", {
#   values <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
#   mock_rast <- create_mock_raster(values)
#   result <- remove_outliers(mock_rast, pmin = 0, pmax = 100)
#   result_values <- values(result)
#
#   # Expect all values to remain NA
#   expect_equal(result_values, values)
# })
#
# test_that("remove_outliers maintains NA values", {
#   values <- c(1, NA, 3, 4, NA, 96, 97, 98, NA, 100)
#   mock_rast <- create_mock_raster(values)
#   result <- remove_outliers(mock_rast)
#   result_values <- values(result)
#
#   # Expect NA values to remain, with outliers set to NA
#   expected <- c(NA, NA, 3, 4, NA, 96, NA, NA, NA, NA)
#   expect_equal(result_values, expected)
# })

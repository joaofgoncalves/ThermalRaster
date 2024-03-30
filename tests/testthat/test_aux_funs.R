
test_that("Check pretty duration output",{

  ti <- Sys.time()
  Sys.sleep(3)
  tf <- Sys.time()
  dt <- tf - ti

  expect_type(pretty_duration(dt), "character")
  expect_equal(pretty_duration(dt), "0 minutes 3 seconds")

})

test_that("returns the last directory name in a standard path", {
  expect_equal(get_sample_dir_name("home/user/documents/project"), "project")
})

test_that("works with a root directory", {
  expect_equal(get_sample_dir_name("/"), "")
})

test_that("works with names without separators", {
  expect_equal(get_sample_dir_name("filename"), "filename")
})

test_that("handles trailing slashes", {
  expect_equal(get_sample_dir_name("home/user/documents/"), "documents")
})

test_that("handles paths with special characters", {
  expect_equal(get_sample_dir_name("home/user/docs/@project!"), "@project!")
  expect_equal(get_sample_dir_name("C:/Program Files/"), "Program Files")
})

test_that("correctly converts 0 to 0", {
  expect_equal(decimal_to_8bit(0), 0)
})

test_that("correctly converts 1 to 255", {
  expect_equal(decimal_to_8bit(1), 255)
})

test_that("correctly converts 0.5 to 128", {
  expect_equal(decimal_to_8bit(0.5), 128)
})

test_that("rounds values correctly", {
  expect_equal(decimal_to_8bit(0.25), 64)
  expect_equal(decimal_to_8bit(0.75), 191)
})

test_that("adds trailing slash if missing", {
  expect_equal(check_trailing_slash("home/user/documents"), "home/user/documents/")
})

test_that("does not add trailing slash if present", {
  expect_equal(check_trailing_slash("home/user/documents/"), "home/user/documents/")
})

test_that("works with root directory", {
  expect_equal(check_trailing_slash("/"), "/")
})

test_that("works with empty string", {
  expect_equal(check_trailing_slash(""), "/")
})

test_that("works with special characters in path", {
  expect_equal(check_trailing_slash("C:/Program Files"), "C:/Program Files/")
})

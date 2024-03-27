
test_that("Check pretty duration output",{

  ti <- Sys.time()
  Sys.sleep(3)
  tf <- Sys.time()
  dt <- tf - ti

  expect_type(pretty_duration(dt), "character")
  expect_equal(pretty_duration(dt), "0 minutes 3 seconds")

})

artifact_path <- function(source_path, installed_path) {
  candidates <- c(
    testthat::test_path("..", "..", source_path),
    file.path(getwd(), source_path),
    system.file(installed_path, package = "eigencore")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) "" else candidates[[1L]]
}

test_that("package timing works with optional memory profiling", {
  skip_if_not_installed("bench")

  timed <- eigencore:::time_repeated(2L, sqrt(4))

  expect_length(timed$times, 2L)
  expect_true(all(is.finite(timed$times)))
  expect_equal(timed$value, 2)
  expect_length(timed$mem_alloc, 1L)
  if (!isTRUE(capabilities("profmem"))) {
    expect_true(is.na(timed$mem_alloc))
  }

  implementation <- paste(deparse(body(eigencore:::time_repeated)), collapse = "\n")
  expect_match(implementation, 'capabilities\\("profmem"\\)')
})

test_that("shipped benchmark entry points do not force memory profiling", {
  artifacts <- c(
    artifact_path(
      file.path("inst", "benchmarks", "_helpers.R"),
      file.path("benchmarks", "_helpers.R")
    ),
    artifact_path(
      file.path("inst", "benchmarks", "bench-readme.R"),
      file.path("benchmarks", "bench-readme.R")
    ),
    artifact_path(
      file.path("vignettes", "benchmarks.Rmd"),
      file.path("doc", "benchmarks.Rmd")
    )
  )

  expect_true(all(nzchar(artifacts)))

  code <- unlist(lapply(artifacts, readLines, warn = FALSE), use.names = FALSE)
  code <- sub("#.*$", "", code)

  expect_false(any(grepl("\\bmemory\\s*=\\s*TRUE\\b", code)))
  expect_gte(
    sum(grepl('capabilities\\s*\\(\\s*["\']profmem["\']', code)),
    length(artifacts)
  )
})

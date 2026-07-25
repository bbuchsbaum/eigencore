local({
  original_capabilities <- base::capabilities
  unlockBinding("capabilities", baseenv())
  assign(
    "capabilities",
    function(what = NULL) {
      if (identical(what, "profmem")) {
        return(FALSE)
      }
      original_capabilities(what)
    },
    envir = baseenv()
  )
  lockBinding("capabilities", baseenv())
  on.exit({
    unlockBinding("capabilities", baseenv())
    assign("capabilities", original_capabilities, envir = baseenv())
    lockBinding("capabilities", baseenv())
  })

  stopifnot(identical(capabilities("profmem"), FALSE))
  cat("Forced capabilities(\"profmem\"):", capabilities("profmem"), "\n")

  library(eigencore)

  package_timing <- eigencore:::time_repeated(2L, sqrt(4))
  stopifnot(
    identical(package_timing$value, 2),
    length(package_timing$times) == 2L,
    all(is.finite(package_timing$times)),
    is.na(package_timing$mem_alloc)
  )

  helper_env <- new.env(parent = globalenv())
  sys.source("inst/benchmarks/_helpers.R", envir = helper_env)
  helper_timing <- helper_env$run_timed(sqrt(4), iterations = 2L)
  stopifnot(
    is.finite(helper_timing$median_seconds),
    is.na(helper_timing$mem_alloc)
  )

  testthat::test_file(
    "tests/testthat/test-profmem-portability.R",
    reporter = "stop"
  )

  output <- rmarkdown::render(
    "vignettes/benchmarks.Rmd",
    output_file = "benchmarks-no-profmem.html",
    output_dir = tempdir(),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  stopifnot(file.exists(output))
  cat("Rendered no-profmem benchmark vignette:", output, "\n")
})

test_that("C++ standard-algorithm users include algorithm directly", {
  src_dir <- testthat::test_path("..", "..", "src")
  sources <- list.files(
    src_dir,
    pattern = "[.](cpp|h|hpp)$",
    full.names = TRUE
  )

  algorithm_call <- paste0(
    "std::(fill|sort|stable_sort|lower_bound|upper_bound|find|find_if|",
    "max_element|min_element|copy|copy_n|transform|reverse|rotate|",
    "partition|nth_element|swap_ranges)[[:space:]]*[(]"
  )
  algorithm_header <-
    "^[[:space:]]*#[[:space:]]*include[[:space:]]*<algorithm>"

  missing_header <- vapply(sources, function(path) {
    lines <- readLines(path, warn = FALSE)
    any(grepl(algorithm_call, lines)) &&
      !any(grepl(algorithm_header, lines))
  }, logical(1))

  expect_equal(
    basename(sources[missing_header]),
    character(),
    info = "Each C++ file that uses a standard algorithm must include <algorithm> directly"
  )
})

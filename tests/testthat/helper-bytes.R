# Compare two raw vectors (e.g. serialize() output) and, on mismatch, report
# only the lengths and the first differing offset instead of a full diff,
# which can run to hundreds of thousands of lines and hide other failures.
expect_same_bytes <- function(actual, expected) {
  same <- identical(actual, expected)
  if (!same) {
    n <- min(length(actual), length(expected))
    diff_at <- which(actual[seq_len(n)] != expected[seq_len(n)])
    first <- if (length(diff_at)) diff_at[[1L]] else n + 1L
    msg <- sprintf(
      "bytes differ: lengths %d vs %d, first difference at offset %d (%d differing bytes in common prefix)",
      length(actual), length(expected), first, length(diff_at)
    )
  } else {
    msg <- ""
  }
  expect(same, msg)
  invisible(actual)
}

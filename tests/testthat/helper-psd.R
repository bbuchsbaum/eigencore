expect_psd_condition <- function(object, class, code, field = NULL) {
  error <- tryCatch(object, error = identity)
  expect_s3_class(error, class)
  expect_s3_class(error, "eigencore_psd_error")
  expect_identical(error$code, code)
  if (!is.null(field)) expect_identical(error$field, field)
  expect_named(
    error,
    c(
      "message", "call", "code", "field", "expected", "actual",
      "source_identity", "factor_identity", "representation", "capability",
      "evidence", "scale", "threshold", "defect", "indices", "details"
    )
  )
  error
}

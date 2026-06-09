#!/usr/bin/env Rscript

source("inst/benchmarks/post-v1-gate-manifest.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

arg_value <- function(args, prefix) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(NULL)
  sub(paste0("^", prefix), "", hit[[1L]])
}

arg_csv <- function(args, prefix) {
  value <- arg_value(args, prefix)
  if (is.null(value) || !nzchar(value)) return(NULL)
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

runner_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  tier <- arg_value(args, "--tier=") %||% "smoke"
  if (!tier %in% c("smoke", "strict", "long")) {
    stop("--tier must be one of smoke, strict, or long.", call. = FALSE)
  }
  list(
    tier = tier,
    gates = arg_csv(args, "--gates="),
    dry_run = "--dry-run" %in% args,
    list_only = "--list" %in% args,
    load_all = "--load-all" %in% args
  )
}

gate_id <- function(gate) gate[["id"]]

select_gate_ids <- function(manifest, tier, selected = NULL) {
  all_ids <- vapply(manifest$gates, gate_id, character(1))
  if (!is.null(selected) && length(selected)) {
    missing <- setdiff(selected, all_ids)
    if (length(missing)) {
      stop(
        "--gates contains unknown gate id(s): ", paste(missing, collapse = ", "),
        ". Available gates: ", paste(all_ids, collapse = ", "),
        call. = FALSE
      )
    }
    return(selected)
  }
  defaults <- manifest$tier_profile[[tier]]$default_gate_ids
  if (identical(defaults, "all")) all_ids else defaults
}

validate_gate_manifest <- function(manifest) {
  required_manifest <- c(
    "version", "generated_on", "artifact_policy", "tier_profile",
    "required_metrics", "current_gate_owner_issue_ids", "gates"
  )
  missing_manifest <- setdiff(required_manifest, names(manifest))
  if (length(missing_manifest)) {
    stop(
      "Post-V1 gate manifest is missing top-level field(s): ",
      paste(missing_manifest, collapse = ", "),
      call. = FALSE
    )
  }
  for (tier in c("smoke", "strict", "long")) {
    if (is.null(manifest$tier_profile[[tier]])) {
      stop("Post-V1 gate manifest is missing tier profile: ", tier, call. = FALSE)
    }
  }

  required_gate <- c(
    "id", "owner_issue", "surface", "script", "command",
    "quick_smoke_command", "cases", "baselines", "artifacts", "thresholds"
  )
  ids <- character()
  for (gate in manifest$gates) {
    missing_gate <- setdiff(required_gate, names(gate))
    if (length(missing_gate)) {
      stop(
        "Gate ", gate$id %||% "<unknown>",
        " is missing field(s): ", paste(missing_gate, collapse = ", "),
        call. = FALSE
      )
    }
    if (gate$id %in% ids) {
      stop("Duplicate post-V1 gate id: ", gate$id, call. = FALSE)
    }
    ids <- c(ids, gate$id)
    if (!gate$owner_issue %in% manifest$current_gate_owner_issue_ids) {
      stop(
        "Gate ", gate$id,
        " owner_issue is not a current gate owner: ", gate$owner_issue,
        call. = FALSE
      )
    }
    if (!file.exists(gate$script)) {
      stop("Gate ", gate$id, " script does not exist: ", gate$script, call. = FALSE)
    }
    if (!length(gate$cases)) {
      stop("Gate ", gate$id, " has no cases.", call. = FALSE)
    }
    if (!length(gate$baselines)) {
      stop("Gate ", gate$id, " has no baselines.", call. = FALSE)
    }
    if (!all(grepl("[.]rds$", gate$artifacts))) {
      stop("Gate ", gate$id, " has non-RDS artifact names.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

gate_command <- function(gate, tier) {
  if (identical(tier, "smoke")) {
    return(gate$quick_smoke_command)
  }
  if (identical(tier, "long")) {
    return(gate$long_command %||% sub("--iterations=3", "--iterations=10", gate$command, fixed = TRUE))
  }
  gate$command
}

command_with_load_all <- function(command) {
  if (!startsWith(command, "Rscript ")) {
    return(command)
  }
  parts <- strsplit(command, " ", fixed = TRUE)[[1L]]
  if (length(parts) < 2L || !grepl("[.]R$", parts[[2L]])) {
    return(command)
  }
  script <- parts[[2L]]
  rest <- if (length(parts) > 2L) paste(parts[-c(1L, 2L)], collapse = " ") else ""
  paste(
    "Rscript -e 'pkgload::load_all(\".\"); source(\"",
    script,
    "\")' --args ",
    rest,
    sep = ""
  )
}

gate_subject <- function(command) {
  hit <- regmatches(command, regexpr("--subject=[^ ]+", command))
  if (!length(hit) || !nzchar(hit)) return(NA_character_)
  sub("^--subject=", "", hit)
}

describe_gate <- function(gate, tier) {
  command <- gate_command(gate, tier)
  data.frame(
    gate_id = gate$id,
    owner_issue = gate$owner_issue,
    tier = tier,
    surface = gate$surface,
    script = gate$script,
    subject = gate_subject(command),
    cases = paste(gate$cases, collapse = ","),
    baselines = paste(gate$baselines, collapse = ","),
    artifacts = paste(gate$artifacts, collapse = ","),
    command = command,
    stringsAsFactors = FALSE
  )
}

run_gate_command <- function(gate, tier, load_all = FALSE) {
  command <- gate_command(gate, tier)
  if (isTRUE(load_all)) {
    command <- command_with_load_all(command)
  }
  message("::group::post-v1 gate ", gate$id)
  message("tier: ", tier)
  message("owner_issue: ", gate$owner_issue)
  message("surface: ", gate$surface)
  message("cases: ", paste(gate$cases, collapse = ", "))
  message("baselines: ", paste(gate$baselines, collapse = ", "))
  message("artifacts: ", paste(gate$artifacts, collapse = ", "))
  message("command: ", command)
  status <- system(command)
  message("::endgroup::")
  if (!identical(status, 0L)) {
    stop(
      "Post-V1 gate failed: ", gate$id,
      " tier=", tier,
      " owner_issue=", gate$owner_issue,
      " exit_status=", status,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

args <- runner_args()
validate_gate_manifest(post_v1_gate_manifest)

selected_ids <- select_gate_ids(post_v1_gate_manifest, args$tier, args$gates)
all_ids <- vapply(post_v1_gate_manifest$gates, gate_id, character(1))
selected <- post_v1_gate_manifest$gates[match(selected_ids, all_ids)]
summary <- do.call(rbind, lapply(selected, describe_gate, tier = args$tier))
row.names(summary) <- NULL
print(summary)

if (args$list_only || args$dry_run) {
  quit(save = "no", status = 0L)
}

dir.create(post_v1_gate_manifest$artifact_policy$directory, recursive = TRUE, showWarnings = FALSE)
for (gate in selected) {
  run_gate_command(gate, args$tier, load_all = args$load_all)
}

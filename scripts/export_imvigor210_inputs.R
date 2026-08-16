#!/usr/bin/env Rscript

# Recreate the two independently verified IMvigor210 inputs used by this project.
# The source is the official IMvigor210CoreBiologies 1.0.0 package archive.
# Patient/sample-level exports are written only to the ignored local data/raw/
# directory and are never prepared for Git or the DOI-backed code archive.

EXPECTED_PACKAGE <- list(
  filename = "IMvigor210CoreBiologies_1.0.0.tar.gz",
  size_bytes = 122127298,
  sha256 = "cfdd3176d7b34de5b04fb9416bfd2b20fa4b6e238aaad5f20b048a34329ea178"
)

EXPECTED_EXPORTS <- list(
  clinical = list(
    filename = "IMvigor210_clinical.csv",
    size_bytes = 76315,
    sha256 = "95b04cb21c885f4b49e1f683d95d44e81a4058e1fd24f3cd15630141292210c0"
  ),
  expression = list(
    filename = "IMvigor210_expression_log2CPM.csv",
    size_bytes = 121488490,
    sha256 = "8cf708b412a4843ad0c2d6d73f3fe804b3308becae8693731113a7b2a573b81c"
  )
)

# The package archive and clinical export retain exact byte gates. Expression
# identity is scientific rather than textual: every in-memory value is checked
# against the direct CPM formula, every value is checked after CSV write/read,
# and the file is checked by the versioned fixed6 semantic contract. The fixed7
# and fixed8 hashes are diagnostics only. No identifier or expression value is
# written to verification reports.
SEMANTIC_CONTRACT_DEFAULT <- file.path(
  "resources", "IMvigor210_expression_semantic_contract_v1.json"
)
SEMANTIC_VERIFIER_DEFAULT <- file.path(
  "scripts", "verify_imvigor210_expression.py"
)
ALL_CELL_ABSOLUTE_TOLERANCE <- 5e-13
ALL_CELL_RELATIVE_TOLERANCE <- 5e-14

# write.table/write.csv use scipen when choosing decimal or scientific form.
# Pin the default used for the canonical exports even if a user profile changes
# it. OutDec covers any classed numeric columns in the clinical data frame.
options(scipen = 0, OutDec = ".")

usage <- function(status = 0L) {
  stream <- if (status == 0L) stdout() else stderr()
  cat(
    paste0(
      "Usage:\n",
      "  Rscript scripts/export_imvigor210_inputs.R \\\n",
      "    --package-tarball /path/to/IMvigor210CoreBiologies_1.0.0.tar.gz \\\n",
      "    [--output-dir data/raw] \\\n",
      "    [--diagnostics-path /path/to/export_diagnostics.tsv] \\\n",
      "    [--semantic-report-path /path/to/expression_semantics.json] \\\n",
      "    [--semantic-contract resources/IMvigor210_expression_semantic_contract_v1.json] \\\n",
      "    [--semantic-verifier scripts/verify_imvigor210_expression.py] \\\n",
      "    [--external-semantic-file-verification]\n",
      "  Rscript scripts/export_imvigor210_inputs.R --verify-only \\\n",
      "    [--output-dir data/raw]\n"
    ),
    file = stream
  )
  quit(save = "no", status = status)
}

parse_args <- function(args) {
  parsed <- list(
    package_tarball = NULL,
    output_dir = file.path("data", "raw"),
    diagnostics_path = NULL,
    semantic_report_path = NULL,
    semantic_contract = SEMANTIC_CONTRACT_DEFAULT,
    semantic_verifier = SEMANTIC_VERIFIER_DEFAULT,
    external_semantic_file_verification = FALSE,
    verify_only = FALSE
  )
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("-h", "--help")) {
      usage(0L)
    } else if (argument == "--verify-only") {
      parsed$verify_only <- TRUE
      index <- index + 1L
    } else if (argument == "--external-semantic-file-verification") {
      parsed$external_semantic_file_verification <- TRUE
      index <- index + 1L
    } else if (argument %in% c(
      "--package-tarball", "--output-dir", "--diagnostics-path",
      "--semantic-report-path", "--semantic-contract", "--semantic-verifier"
    )) {
      if (index == length(args)) {
        stop("Missing value after ", argument, call. = FALSE)
      }
      value <- args[[index + 1L]]
      if (argument == "--package-tarball") {
        parsed$package_tarball <- value
      } else if (argument == "--diagnostics-path") {
        parsed$diagnostics_path <- value
      } else if (argument == "--semantic-report-path") {
        parsed$semantic_report_path <- value
      } else if (argument == "--semantic-contract") {
        parsed$semantic_contract <- value
      } else if (argument == "--semantic-verifier") {
        parsed$semantic_verifier <- value
      } else {
        parsed$output_dir <- value
      }
      index <- index + 2L
    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
  }
  if (!parsed$verify_only && is.null(parsed$package_tarball)) {
    stop("--package-tarball is required unless --verify-only is used", call. = FALSE)
  }
  if (parsed$verify_only && parsed$external_semantic_file_verification) {
    stop(
      "--external-semantic-file-verification cannot be used with --verify-only",
      call. = FALSE
    )
  }
  parsed
}

sha256_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(tolower(digest::digest(
      object = path, algo = "sha256", serialize = FALSE, file = TRUE
    )))
  }

  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    output <- system2(
      sha256sum, args = shQuote(path), stdout = TRUE, stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status) || identical(status, 0L)) {
      digest <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
      if (grepl("^[0-9A-Fa-f]{64}$", digest)) {
        return(tolower(digest))
      }
    }
  }

  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    output <- system2(
      shasum, args = c("-a", "256", shQuote(path)),
      stdout = TRUE, stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status) || identical(status, 0L)) {
      digest <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
      if (grepl("^[0-9A-Fa-f]{64}$", digest)) {
        return(tolower(digest))
      }
    }
  }

  stop(
    "SHA-256 verification requires the R package 'digest', sha256sum, or shasum.",
    call. = FALSE
  )
}

compare_expression_matrices <- function(observed, expected, label) {
  if (!identical(dim(observed), dim(expected))) {
    stop(label, " dimensions differ.", call. = FALSE)
  }
  if (!identical(dimnames(observed), dimnames(expected))) {
    stop(label, " ordered feature or sample identifiers differ.", call. = FALSE)
  }
  if (
    anyNA(observed) || anyNA(expected) ||
      any(!is.finite(observed)) || any(!is.finite(expected))
  ) {
    stop(label, " contains a non-finite value.", call. = FALSE)
  }
  if (any(observed < 0) || any(expected < 0)) {
    stop(label, " contains a negative value.", call. = FALSE)
  }
  if (!identical(observed == 0, expected == 0)) {
    stop(label, " structural-zero masks differ.", call. = FALSE)
  }

  absolute_delta <- abs(observed - expected)
  tolerance <- ALL_CELL_ABSOLUTE_TOLERANCE +
    ALL_CELL_RELATIVE_TOLERANCE * abs(expected)
  within_tolerance <- all(absolute_delta <= tolerance)
  nonzero <- expected != 0
  max_relative_delta <- if (any(nonzero)) {
    max(absolute_delta[nonzero] / abs(expected[nonzero]))
  } else {
    0
  }
  diagnostics <- list(
    n_rows = nrow(observed),
    n_columns = ncol(observed),
    max_absolute_delta = max(absolute_delta),
    max_relative_delta = max_relative_delta,
    within_tolerance = within_tolerance
  )
  rm(absolute_delta, tolerance, nonzero)
  if (!within_tolerance) {
    stop(label, " exceeds the all-cell numeric tolerance.", call. = FALSE)
  }
  diagnostics
}

read_expression_csv <- function(path) {
  value <- utils::read.csv(
    path,
    header = TRUE,
    row.names = 1L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  value
}

python_executable <- function() {
  candidates <- Sys.which(c("python3", "python"))
  candidates <- unname(candidates[nzchar(candidates)])
  if (length(candidates) == 0L) {
    stop(
      "Python is required for IMvigor210 expression semantic verification.",
      call. = FALSE
    )
  }
  candidates[[1L]]
}

verify_expression_file_semantically <- function(
  path, semantic_contract, semantic_verifier, report_path = NULL
) {
  if (!file.exists(semantic_contract) || dir.exists(semantic_contract)) {
    stop("Semantic contract is missing.", call. = FALSE)
  }
  if (!file.exists(semantic_verifier) || dir.exists(semantic_verifier)) {
    stop("Semantic verifier is missing.", call. = FALSE)
  }
  command_arguments <- c(
    semantic_verifier,
    "--input", path,
    "--contract", semantic_contract
  )
  if (!is.null(report_path)) {
    command_arguments <- c(command_arguments, "--report", report_path)
  }
  output <- system2(
    python_executable(),
    args = shQuote(command_arguments),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      "Expression semantic verifier failed. Its output is redacted; ",
      paste(output, collapse = "\n"),
      call. = FALSE
    )
  }
  message("VERIFIED\timvigor210_expression_semantic_contract_v1")
  invisible(TRUE)
}

write_export_diagnostics <- function(
  path, staged_paths, direct_diagnostics, readback_diagnostics
) {
  if (is.null(path)) {
    return(invisible(NULL))
  }
  if (dir.exists(path)) {
    stop(
      "--diagnostics-path must name a file, not a directory: ", path,
      call. = FALSE
    )
  }
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    stop("Could not create diagnostics directory: ", parent, call. = FALSE)
  }

  clinical_size <- unname(file.info(staged_paths[["clinical"]])$size)
  expression_size <- unname(file.info(staged_paths[["expression"]])$size)
  observed <- c(
    clinical_size_bytes = as.character(clinical_size),
    clinical_sha256 = sha256_file(staged_paths[["clinical"]]),
    expression_size_bytes = as.character(expression_size),
    expression_sha256 = sha256_file(staged_paths[["expression"]]),
    expression_rows = as.character(direct_diagnostics$n_rows),
    expression_columns = as.character(direct_diagnostics$n_columns),
    direct_formula_max_absolute_delta = format(
      direct_diagnostics$max_absolute_delta, scientific = TRUE, digits = 17
    ),
    direct_formula_max_relative_delta = format(
      direct_diagnostics$max_relative_delta, scientific = TRUE, digits = 17
    ),
    direct_formula_all_cells_within_tolerance =
      as.character(direct_diagnostics$within_tolerance),
    write_readback_max_absolute_delta = format(
      readback_diagnostics$max_absolute_delta, scientific = TRUE, digits = 17
    ),
    write_readback_max_relative_delta = format(
      readback_diagnostics$max_relative_delta, scientific = TRUE, digits = 17
    ),
    write_readback_all_cells_within_tolerance =
      as.character(readback_diagnostics$within_tolerance)
  )
  expected <- c(
    clinical_size_bytes = as.character(EXPECTED_EXPORTS$clinical$size_bytes),
    clinical_sha256 = EXPECTED_EXPORTS$clinical$sha256,
    expression_size_bytes = as.character(EXPECTED_EXPORTS$expression$size_bytes),
    expression_sha256 = EXPECTED_EXPORTS$expression$sha256,
    expression_rows = "31286",
    expression_columns = "348",
    direct_formula_max_absolute_delta = "aggregate diagnostic",
    direct_formula_max_relative_delta = "aggregate diagnostic",
    direct_formula_all_cells_within_tolerance = "TRUE",
    write_readback_max_absolute_delta = "aggregate diagnostic",
    write_readback_max_relative_delta = "aggregate diagnostic",
    write_readback_all_cells_within_tolerance = "TRUE"
  )
  required <- c(
    clinical_size_bytes = TRUE,
    clinical_sha256 = TRUE,
    expression_size_bytes = FALSE,
    expression_sha256 = FALSE,
    expression_rows = TRUE,
    expression_columns = TRUE,
    direct_formula_max_absolute_delta = FALSE,
    direct_formula_max_relative_delta = FALSE,
    direct_formula_all_cells_within_tolerance = TRUE,
    write_readback_max_absolute_delta = FALSE,
    write_readback_max_relative_delta = FALSE,
    write_readback_all_cells_within_tolerance = TRUE
  )
  match <- observed == expected
  match[grepl("_max_(absolute|relative)_delta$", names(match))] <- NA
  diagnostics <- data.frame(
    metric = names(observed),
    observed = unname(observed),
    expected = unname(expected),
    required = unname(required),
    match = unname(match),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    diagnostics, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  message("DIAGNOSTICS\t", normalizePath(path, mustWork = TRUE))
  invisible(diagnostics)
}

verify_file <- function(path, specification, label) {
  if (!file.exists(path) || dir.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  observed_size <- unname(file.info(path)$size)
  observed_sha256 <- sha256_file(path)
  message(
    "OBSERVED\t", label, "\tsize_bytes=", observed_size,
    "\tsha256=", observed_sha256
  )
  if (is.na(observed_size) || observed_size != specification$size_bytes) {
    stop(
      label, " size mismatch: ", observed_size, " != ",
      specification$size_bytes, "; observed SHA-256: ", observed_sha256,
      call. = FALSE
    )
  }
  if (!identical(observed_sha256, specification$sha256)) {
    stop(
      label, " SHA-256 mismatch: ", observed_sha256, " != ",
      specification$sha256, call. = FALSE
    )
  }
  message("VERIFIED\t", label, "\t", normalizePath(path, mustWork = TRUE))
  invisible(TRUE)
}

verify_exports <- function(output_dir, arguments, report_path = NULL) {
  verify_file(
    file.path(output_dir, EXPECTED_EXPORTS$clinical$filename),
    EXPECTED_EXPORTS$clinical,
    "imvigor210_clinical_export"
  )
  verify_expression_file_semantically(
    file.path(output_dir, EXPECTED_EXPORTS$expression$filename),
    arguments$semantic_contract,
    arguments$semantic_verifier,
    report_path
  )
  invisible(TRUE)
}

write_canonical_csv <- function(value, path) {
  connection <- file(path, open = "wb")
  tryCatch(
    utils::write.csv(
      value, connection, row.names = TRUE, quote = TRUE,
      na = "NA", eol = "\n"
    ),
    finally = close(connection)
  )
}

publish_verified_exports <- function(staged_paths, output_dir, arguments) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }

  final_paths <- vapply(
    EXPECTED_EXPORTS,
    function(specification) file.path(output_dir, specification$filename),
    character(1L)
  )
  if (file.exists(final_paths[["clinical"]])) {
    verify_file(
      final_paths[["clinical"]], EXPECTED_EXPORTS$clinical, "existing_clinical"
    )
  }
  if (file.exists(final_paths[["expression"]])) {
    verify_expression_file_semantically(
      final_paths[["expression"]],
      arguments$semantic_contract,
      arguments$semantic_verifier,
      arguments$semantic_report_path
    )
  }

  temporary_paths <- stats::setNames(
    paste0(unname(final_paths), ".tmp"),
    names(final_paths)
  )
  on.exit(unlink(temporary_paths[file.exists(temporary_paths)]), add = TRUE)
  for (name in names(final_paths)) {
    if (file.exists(final_paths[[name]])) {
      next
    }
    if (file.exists(temporary_paths[[name]])) {
      unlink(temporary_paths[[name]])
    }
    copied <- file.copy(
      staged_paths[[name]], temporary_paths[[name]], overwrite = FALSE,
      copy.mode = TRUE, copy.date = FALSE
    )
    if (!copied) {
      stop("Could not stage output: ", final_paths[[name]], call. = FALSE)
    }
    if (name == "clinical") {
      verify_file(
        temporary_paths[[name]], EXPECTED_EXPORTS[[name]], "staged_clinical"
      )
    } else {
      if (!identical(
        sha256_file(temporary_paths[[name]]), sha256_file(staged_paths[[name]])
      )) {
        stop("Expression export changed during local publication.", call. = FALSE)
      }
    }
  }

  for (name in names(final_paths)) {
    if (file.exists(final_paths[[name]])) {
      next
    }
    if (!file.rename(temporary_paths[[name]], final_paths[[name]])) {
      stop("Could not publish output: ", final_paths[[name]], call. = FALSE)
    }
  }
  verify_file(
    final_paths[["clinical"]], EXPECTED_EXPORTS$clinical, "published_clinical"
  )
  if (!file.exists(final_paths[["expression"]])) {
    stop("Published expression export is missing.", call. = FALSE)
  }
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
if (arguments$verify_only) {
  verify_exports(
    arguments$output_dir, arguments, arguments$semantic_report_path
  )
  quit(save = "no", status = 0L)
}

package_tarball <- arguments$package_tarball
if (basename(package_tarball) != EXPECTED_PACKAGE$filename) {
  stop(
    "The package archive must retain its canonical filename: ",
    EXPECTED_PACKAGE$filename, call. = FALSE
  )
}
verify_file(package_tarball, EXPECTED_PACKAGE, "imvigor210_processed_package")

required_packages <- c("DESeq", "Biobase", "edgeR")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L)
)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R/Bioconductor package(s): ", paste(missing_packages, collapse = ", "),
    ". Use an R/Bioconductor environment compatible with the legacy DESeq ",
    "CountDataSet stored in IMvigor210CoreBiologies 1.0.0.",
    call. = FALSE
  )
}

staging_dir <- tempfile(pattern = "imvigor210-export-")
dir.create(staging_dir)
on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
archive_dir <- file.path(staging_dir, "archive")
dir.create(archive_dir)
utils::untar(package_tarball, exdir = archive_dir)

data_candidates <- list.files(
  archive_dir,
  pattern = "^cds[.](rda|RData|rdata)$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(data_candidates) != 1L) {
  stop(
    "Expected exactly one cds data object in the verified package; found ",
    length(data_candidates), call. = FALSE
  )
}

object_environment <- new.env(parent = emptyenv())
loaded_names <- load(data_candidates[[1L]], envir = object_environment)
if (!("cds" %in% loaded_names)) {
  stop("The verified package data file did not load an object named 'cds'.", call. = FALSE)
}
cds <- object_environment$cds
count_matrix <- as.matrix(DESeq::counts(cds))
clinical <- as.data.frame(Biobase::pData(cds), stringsAsFactors = FALSE)

if (!identical(dim(count_matrix), c(31286L, 348L))) {
  stop(
    "Unexpected count-matrix dimensions: ", paste(dim(count_matrix), collapse = " x "),
    call. = FALSE
  )
}
if (!identical(dim(clinical), c(348L, 25L))) {
  stop(
    "Unexpected clinical-table dimensions: ", paste(dim(clinical), collapse = " x "),
    call. = FALSE
  )
}
if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix))) {
  stop("The count matrix lacks feature or sample identifiers.", call. = FALSE)
}
if (anyDuplicated(rownames(count_matrix)) || anyDuplicated(colnames(count_matrix))) {
  stop("The count matrix contains duplicate feature or sample identifiers.", call. = FALSE)
}
if (!identical(colnames(count_matrix), rownames(clinical))) {
  stop("Clinical rows and expression columns differ in identity or order.", call. = FALSE)
}
if (anyNA(count_matrix) || any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
  stop("The count matrix contains missing, non-finite or negative values.", call. = FALSE)
}
if (any(count_matrix != floor(count_matrix))) {
  stop("The source object is not an integer count matrix.", call. = FALSE)
}

# This is the historical project transform. edgeR::cpm() is evaluated on the
# complete 31,286-feature count matrix, then a pseudocount of one is added on
# the CPM scale before the base-2 logarithm. No sample or feature is filtered.
expression_log2cpm <- log2(edgeR::cpm(count_matrix, log = FALSE) + 1)
if (!identical(dimnames(expression_log2cpm), dimnames(count_matrix))) {
  stop("Expression normalization changed source identifiers or order.", call. = FALSE)
}
if (anyNA(expression_log2cpm) || any(!is.finite(expression_log2cpm))) {
  stop("The normalized expression matrix contains non-finite values.", call. = FALSE)
}

library_sizes <- colSums(count_matrix)
if (anyNA(library_sizes) || any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
  stop("The count matrix contains an invalid library size.", call. = FALSE)
}
direct_expression_log2cpm <- log2(
  sweep(count_matrix, 2L, library_sizes, "/") * 1e6 + 1
)
direct_diagnostics <- compare_expression_matrices(
  expression_log2cpm,
  direct_expression_log2cpm,
  "edgeR CPM and direct-formula expression"
)
rm(direct_expression_log2cpm, library_sizes)

staged_paths <- c(
  clinical = file.path(staging_dir, EXPECTED_EXPORTS$clinical$filename),
  expression = file.path(staging_dir, EXPECTED_EXPORTS$expression$filename)
)
write_canonical_csv(clinical, staged_paths[["clinical"]])
write_canonical_csv(expression_log2cpm, staged_paths[["expression"]])

readback_expression_log2cpm <- read_expression_csv(staged_paths[["expression"]])
readback_diagnostics <- compare_expression_matrices(
  readback_expression_log2cpm,
  expression_log2cpm,
  "Written/read-back expression and in-memory expression"
)
rm(readback_expression_log2cpm)
write_export_diagnostics(
  arguments$diagnostics_path,
  staged_paths,
  direct_diagnostics,
  readback_diagnostics
)
verify_file(staged_paths[["clinical"]], EXPECTED_EXPORTS$clinical, "generated_clinical")
if (!arguments$external_semantic_file_verification) {
  verify_expression_file_semantically(
    staged_paths[["expression"]],
    arguments$semantic_contract,
    arguments$semantic_verifier,
    arguments$semantic_report_path
  )
} else {
  message(
    "DEFERRED\timvigor210_expression_semantic_contract_v1\t",
    "external caller must verify before use"
  )
}
publish_verified_exports(staged_paths, arguments$output_dir, arguments)

message("R_VERSION\t", R.version.string)
message("DESEQ_VERSION\t", as.character(utils::packageVersion("DESeq")))
message("EDGER_VERSION\t", as.character(utils::packageVersion("edgeR")))

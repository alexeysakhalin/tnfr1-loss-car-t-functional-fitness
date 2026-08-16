#!/usr/bin/env Rscript

# Recreate the two checksum-pinned IMvigor210 inputs used by this project.
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

# Byte-for-byte CSV identity remains the release gate. This second contract is
# deliberately insensitive to harmless numeric text formatting: it hashes the
# complete row-major expression matrix after fixed eight-decimal formatting,
# while feature and sample identifiers are hashed separately as UTF-8 lines.
# The expected values were calculated once from the checksum-locked project
# export. No identifier or expression value is written to diagnostics.
EXPECTED_EXPRESSION_SEMANTICS <- list(
  n_rows = 31286L,
  n_columns = 348L,
  value_format = "%.8f",
  row_ids_sha256 = "99c0a222c27bae6c35d479da88aeb8812d5721875f3d5705c5711bffc48c364d",
  sample_ids_sha256 = "388ef1e09720f61bce1939b15a3b39eec989d415ad590fa07a03dc1ec68619e8",
  values_sha256 = "da0c1007d1a267f86cb4dbfc0dca85eb8204ef7f2439e60b74ef56eea9d92444"
)

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
      "    [--diagnostics-path /path/to/export_diagnostics.tsv]\n",
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
    } else if (argument %in% c(
      "--package-tarball", "--output-dir", "--diagnostics-path"
    )) {
      if (index == length(args)) {
        stop("Missing value after ", argument, call. = FALSE)
      }
      value <- args[[index + 1L]]
      if (argument == "--package-tarball") {
        parsed$package_tarball <- value
      } else if (argument == "--diagnostics-path") {
        parsed$diagnostics_path <- value
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

sha256_utf8_lines <- function(values) {
  if (anyNA(values) || any(grepl("[\r\n]", values))) {
    stop(
      "Identifiers for semantic hashing must be non-missing single lines.",
      call. = FALSE
    )
  }
  path <- tempfile(pattern = "imvigor210-identifiers-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  tryCatch(
    for (value in values) {
      writeBin(charToRaw(enc2utf8(paste0(value, "\n"))), connection)
    },
    finally = close(connection)
  )
  sha256_file(path)
}

expression_semantic_diagnostics <- function(value) {
  if (is.null(rownames(value)) || is.null(colnames(value))) {
    stop(
      "Semantic hashing requires feature and sample identifiers.",
      call. = FALSE
    )
  }
  if (anyNA(value) || any(!is.finite(value))) {
    stop(
      "Semantic hashing requires a finite expression matrix.",
      call. = FALSE
    )
  }

  path <- tempfile(pattern = "imvigor210-values-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  old_numeric_locale <- Sys.getlocale("LC_NUMERIC")
  on.exit(
    suppressWarnings(Sys.setlocale("LC_NUMERIC", old_numeric_locale)),
    add = TRUE
  )
  if (is.na(suppressWarnings(Sys.setlocale("LC_NUMERIC", "C")))) {
    close(connection)
    stop(
      "Could not set the C numeric locale for semantic hashing.",
      call. = FALSE
    )
  }
  tryCatch(
    for (row_index in seq_len(nrow(value))) {
      canonical_row <- paste0(
        paste(
          sprintf(
            EXPECTED_EXPRESSION_SEMANTICS$value_format,
            value[row_index, , drop = TRUE]
          ),
          collapse = ","
        ),
        "\n"
      )
      writeBin(charToRaw(canonical_row), connection)
    },
    finally = close(connection)
  )

  list(
    n_rows = nrow(value),
    n_columns = ncol(value),
    row_ids_sha256 = sha256_utf8_lines(rownames(value)),
    sample_ids_sha256 = sha256_utf8_lines(colnames(value)),
    values_sha256 = sha256_file(path)
  )
}

verify_expression_semantics <- function(observed) {
  expected <- c(
    n_rows = as.character(EXPECTED_EXPRESSION_SEMANTICS$n_rows),
    n_columns = as.character(EXPECTED_EXPRESSION_SEMANTICS$n_columns),
    row_ids_sha256 = EXPECTED_EXPRESSION_SEMANTICS$row_ids_sha256,
    sample_ids_sha256 = EXPECTED_EXPRESSION_SEMANTICS$sample_ids_sha256,
    values_sha256 = EXPECTED_EXPRESSION_SEMANTICS$values_sha256
  )
  observed_vector <- c(
    n_rows = as.character(observed$n_rows),
    n_columns = as.character(observed$n_columns),
    row_ids_sha256 = observed$row_ids_sha256,
    sample_ids_sha256 = observed$sample_ids_sha256,
    values_sha256 = observed$values_sha256
  )
  mismatched <- names(expected)[observed_vector != expected]
  if (length(mismatched) > 0L) {
    stop(
      "Generated expression semantic contract mismatch: ",
      paste(mismatched, collapse = ", "),
      call. = FALSE
    )
  }
  message(
    "VERIFIED\texpression_semantics\tvalues_sha256=",
    observed$values_sha256
  )
  invisible(TRUE)
}

write_export_diagnostics <- function(path, staged_paths, expression_semantics) {
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
    expression_rows = as.character(expression_semantics$n_rows),
    expression_columns = as.character(expression_semantics$n_columns),
    expression_row_ids_sha256 = expression_semantics$row_ids_sha256,
    expression_sample_ids_sha256 = expression_semantics$sample_ids_sha256,
    expression_values_8dp_sha256 = expression_semantics$values_sha256
  )
  expected <- c(
    clinical_size_bytes = as.character(EXPECTED_EXPORTS$clinical$size_bytes),
    clinical_sha256 = EXPECTED_EXPORTS$clinical$sha256,
    expression_size_bytes = as.character(EXPECTED_EXPORTS$expression$size_bytes),
    expression_sha256 = EXPECTED_EXPORTS$expression$sha256,
    expression_rows = as.character(EXPECTED_EXPRESSION_SEMANTICS$n_rows),
    expression_columns = as.character(EXPECTED_EXPRESSION_SEMANTICS$n_columns),
    expression_row_ids_sha256 = EXPECTED_EXPRESSION_SEMANTICS$row_ids_sha256,
    expression_sample_ids_sha256 =
      EXPECTED_EXPRESSION_SEMANTICS$sample_ids_sha256,
    expression_values_8dp_sha256 =
      EXPECTED_EXPRESSION_SEMANTICS$values_sha256
  )
  diagnostics <- data.frame(
    metric = names(observed),
    observed = unname(observed),
    expected = unname(expected),
    match = unname(observed == expected),
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

verify_exports <- function(output_dir) {
  for (name in names(EXPECTED_EXPORTS)) {
    specification <- EXPECTED_EXPORTS[[name]]
    verify_file(
      file.path(output_dir, specification$filename),
      specification,
      paste0("imvigor210_", name, "_export")
    )
  }
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

publish_verified_exports <- function(staged_paths, output_dir) {
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
  for (name in names(final_paths)) {
    final_path <- final_paths[[name]]
    if (file.exists(final_path)) {
      verify_file(final_path, EXPECTED_EXPORTS[[name]], paste0("existing_", name))
    }
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
    verify_file(
      temporary_paths[[name]], EXPECTED_EXPORTS[[name]], paste0("staged_", name)
    )
  }

  for (name in names(final_paths)) {
    if (file.exists(final_paths[[name]])) {
      next
    }
    if (!file.rename(temporary_paths[[name]], final_paths[[name]])) {
      stop("Could not publish output: ", final_paths[[name]], call. = FALSE)
    }
  }
  verify_exports(output_dir)
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
if (arguments$verify_only) {
  verify_exports(arguments$output_dir)
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

staged_paths <- c(
  clinical = file.path(staging_dir, EXPECTED_EXPORTS$clinical$filename),
  expression = file.path(staging_dir, EXPECTED_EXPORTS$expression$filename)
)
write_canonical_csv(clinical, staged_paths[["clinical"]])
write_canonical_csv(expression_log2cpm, staged_paths[["expression"]])

expression_semantics <- expression_semantic_diagnostics(expression_log2cpm)
write_export_diagnostics(
  arguments$diagnostics_path, staged_paths, expression_semantics
)
verify_expression_semantics(expression_semantics)
verify_file(staged_paths[["clinical"]], EXPECTED_EXPORTS$clinical, "generated_clinical")
verify_file(
  staged_paths[["expression"]], EXPECTED_EXPORTS$expression,
  "generated_expression"
)
publish_verified_exports(staged_paths, arguments$output_dir)

message("R_VERSION\t", R.version.string)
message("DESEQ_VERSION\t", as.character(utils::packageVersion("DESeq")))
message("EDGER_VERSION\t", as.character(utils::packageVersion("edgeR")))

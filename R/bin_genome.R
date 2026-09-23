.binGenome_env <- new.env(parent = emptyenv())
.binGenome_env$binGenome_path <- NULL # NULL = use the genomeBinner.jar shipped with the package

#' Wrapper for binGenome module
#'
#' @rdname bin_genome
#' @name bin_genome
#'
#' @param path Path to a binGenome.sh executable or a genomeBinner jar file.
#'  Only needed to override the default, the genomeBinner.jar shipped with the package.
#'
#' @return Invisibly the path on successful execution.
#' @export
set_binGenome_path <- function(path) {
  if (!file.exists(path)) stop("Invalid binGenome.sh path: ", path)
  if (!endsWith(path, ".sh") && !endsWith(path, ".jar")) stop("The file should be a shell executable (.sh) or a jar file (.jar).")
  .binGenome_env$binGenome_path <- normalizePath(path)
  invisible(path)
}

#' @rdname bin_genome
#' @description
#' Runs the binGenome.sh.
#' Prints the log of how many regions were too small for binning in each sample.
#' If stranded, pos and neg are merged and and named after the file given as positive strand (minus/neg for strand -1)
#'
#' @param input_dir A character string of path to the folder containing bedgraph files.
#' @param strand An integer in c(-1, 0, 1) to specify if and how the strandedness
#' of coverage files should be considered.
#' 0: no filter | 1: regular stranded | -1: reverse stranded.
#' if `strand != 0`, it groups files by applying grep according to
#' the `strand_pattern`.
#' @param annotation A character string of path to annotation file in BED format.
#' @param output_dir Output folder for binGenome.sh
#' @param strand_pattern Identifiers for positive and negative strand to match
#' in the names of annotation files, separated by '/'.
#' First element corresponds to positive and second to negative. Default is "pos/neg".
#' @param bins Number of bins to bin the gene body.
#' @param grep_pattern Filter for a grep pattern in files. perl is TRUE.
#' @param fixedBinSizeUpstream,fixedBinSizeDownstream Specifies the parameters that will be
#' used for binning the downstream or upstream regions, respectively. Format: "binsize:binnumber".
#' @param bedgraphNames,annotationNames A character vector of names for bedgraph and annotation files.
#' @param cores Number of threads to use. Default: max(1, parallel::detectCores() - 2)
#' @param normalize Logical. If `TRUE`, an additional per-gene normalized version is written.
#' @param tmpDir Path to a tmp directory. If NULL (default), java temp dir is used. Set a tempDir if the default is too small or not writable
#'
#' @returns Returns the `output_dir` invisibly.
#'
#' @details
#' For binning, only files with .bedgraph or .bedgraph.gz extension in the `input_dir` are considered.
#' The function filters the files by "\\.bedgraph$", hence grep_pattern does not have to include the extension.
#'
#' The `strand_pattern` is used literally, so user should not escape regex metacharacters.
#' For example, by strand_pattern = "+/-", plus is matched literally without any problem.
#'
#' Jar file's output is saved to the `out/` directory inside `output_dir`.
#'
#' Since user quotas might be restricted on clusters, the jar file might often
#' crash due to low space. This is why it might help to set `-Djava.io.tmpdir`
#' to a directory where there is more allowance.
#'
#' @export
bin_genome <- function(input_dir, strand, annotation, output_dir,
                       strand_pattern = "pos/neg",
                       bins = NULL, grep_pattern = NULL,
                       fixedBinSizeUpstream = NULL, fixedBinSizeDownstream = NULL,
                       bedgraphNames = NULL, annotationNames = NULL,
                       cores = NULL, normalize = FALSE, tmpDir = NULL) {
  if (is.null(.binGenome_env$binGenome_path)) {
    .binGenome_env$binGenome_path <- system.file("binGenome_environment", "genomeBinner.jar", package = "metageneR", mustWork = TRUE)
  }
  message(paste("Using programme:", .binGenome_env$binGenome_path))

  if (endsWith(.binGenome_env$binGenome_path, ".jar") && is.null(tmpDir)) {
    message("No tmpDir set, the Java default temp directory is used. Set tmpDir if it is small or not writable.")
  }

  # Argument checks
  if (!dir.exists(input_dir)) {
    stop("input_dir does not exist: ", input_dir)
  }
  if (!file.exists(annotation)) {
    stop("annotation file does not exist: ", annotation)
  }
  if (is.null(bins)) {
    stop("You must provide `bins`")
  }
  # 0: not stranded, 1: regular, -1: reverse stranded
  if (!(strand %in% c(-1, 0, 1))) {
    stop(paste("Unrecognized argument for strand:", strand,
               "\nShould be one of c(-1, 0, 1)"))
  }
  if (!grepl("^[^/]+/[^/]+$", strand_pattern)) {
    stop(paste("Wrong format for strand_pattern argument: ", strand_pattern,
               "\nShould be two words separated by exactly one '/'."))
  }

  # Create the output directory if it doesn't exist
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  } else if (dir.exists(output_dir)) {
    if (interactive()) {
      should_overwrite <- askYesNo("Output directory exists. Should its content be overwritten?",
                                   prompts = getOption("askYesNo", gettext(c("yes", "no", "Cancel"))))
      if (is.na(should_overwrite) || !should_overwrite) stop("Aborting...")
    } else {
      warning("Overwriting the output directory: ", output_dir)
    }
  }

  # if cores not given, automatically detect
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 2)
  }
  params <- as.list(environment())

  # =========== Config and Logs  ===========
  write_config("bin_genome", path = file.path(output_dir, "config.json"), params)
  # create logger function
  logger <- create_logger("bin_genome", params, log_file = file.path(output_dir, "log.txt"))
  on.exit(logger(close = TRUE), add = T)
  # ================================

  all_files <- list.files(input_dir)
  if (!is.null(grep_pattern)) {
    all_files <- grep(grep_pattern, all_files, perl = TRUE, value = TRUE)
  }

  # Process files: for each .gz file, check if uncompressed counterpart exists
  files <- character()
  for (f in all_files) {
    if (grepl("\\.bedgraph\\.gz$", f, ignore.case = TRUE)) {
      base_name <- sub("\\.gz$", "", f, ignore.case = TRUE)
      if (file.exists(file.path(input_dir, base_name))) {
        files <- c(files, base_name)
      } else {
        stop("Found .bedgraph.gz file '", f, "' but its uncompressed counterpart '", base_name,
             "' does not exist. Please decompress it first using decompress_bedgraph_gz() function.")
      }
    } else if (grepl("\\.bedgraph$", f, ignore.case = TRUE)) {
      files <- c(files, f)
    }
  }

  if (length(files) == 0) {
    stop("No .bedgraph files found in ", input_dir)
  }

  # ensure chr_name/seq_name in BED and BEDGRAPH is the same
  validate_seqnames(annotation, file.path(input_dir, files), logger)

  # If stranded, group files
  pos <- neg <- bedgraph <- bedgraphPos <- bedgraphNeg <- NULL
  strand_words <- strsplit(strand_pattern, "/", fixed = TRUE)[[1]]
  if (strand != 0) {
    if (strand == -1) {
      pos <- strand_words[2]
      neg <- strand_words[1]
    } else if (strand == 1) {
      pos <- strand_words[1]
      neg <- strand_words[2]
    }
    bedgraphPos <- paste(file.path(input_dir, grep(pos, files, value = T, fixed = T)), collapse = ",")
    bedgraphNeg <- paste(file.path(input_dir, grep(neg, files, value = T, fixed = T)), collapse = ",")
  } else {
    bedgraph <- paste(file.path(input_dir, files), collapse = ",")
  }


  # save arguments in a named list, NULL is okay
  opts <- list(
    "--annotation" = annotation,
    "--bedgraph" = bedgraph,
    "--bedgraphNeg" = bedgraphNeg,
    "--bedgraphPos" = bedgraphPos,
    "--outputDir" = file.path(output_dir, "out"),
    "--bins" = bins,
    "--fixedBinSizeUpstream" = fixedBinSizeUpstream,
    "--fixedBinSizeDownstream" = fixedBinSizeDownstream,
    "--bedgraphNames" = bedgraphNames,
    "--annotationNames" = annotationNames,
    "--cores" = cores,
    "--tmpDir" = tmpDir
    )

  # flatten -> c("--flag", "val", ...) and remove the NULL args
  args <- unlist(Map(
    function(flag, val) {
      if (!is.null(val)) c(flag, as.character(val)) else NULL
      }, names(opts), opts), use.names = FALSE)
  # add normalize if given
  if (normalize) {
    args <- c(args, "--normalize")
  }

  executable <- .binGenome_env$binGenome_path
  # Handle jar file case: 1) change executable 2) specify -Djava.io.tmpdir
  if (endsWith(executable, ".jar")) {
    executable <- "java"
    args <- c("-jar", .binGenome_env$binGenome_path, args)
    if (!is.null(tmpDir)) args <- c(paste0("-Djava.io.tmpdir=", tmpDir), args)
  }

  # Run process
  console_output <- processx::run(
    command = executable,
    args = args,
    echo = FALSE,
    echo_cmd = TRUE,
    stdout_line_callback = function(line, proc) logger(line, type = "message"),
    stderr_line_callback = function(line, proc) logger(line, type = "warning")
    # line callbacks write clutters the line beginnings with [warning] but is fine...
    # below code will not work due to some weird error, maybe cos stderr too long?..
  )

  # if (console_output$status != 0) {
  #   message <- paste(sprintf("processx::run failed with status %d", console_output$status),
  #                    "Stderr:", console_output$stderr,
  #                    "Stdout:", console_output$stdout, collapse = "\n")
  #   logger(message, type = "error")
  # }

  warns <- too_small_warns_by_sample(console_output$stdout, strand)
  if (length(warns) != 0) {
    logger("Too small regions found:")
    for (sample in names(warns)) {
      logger(sprintf("  %s: %d", sample, warns[[sample]]))
    }
  }

  invisible(output_dir)
}

### HELPERS ===============

too_small_warns_by_sample <- function(console_output, strand) {
  sample <- NULL
  log_by_sample <- list() # list will be named with samples/files
  log_lines <- strsplit(console_output, split = "\n")[[1]]

  for (line in log_lines) {
    # Start of a new file
    if (grepl("^\\[INFO\\] started with.*\\.bedgraph'", line)) {
      # Extract filename from path
      path_match <- regmatches(line, regexpr("'[^']+\\.bedgraph'", line))
      sample <- gsub(".*\\/|\\.bedgraph'", "", path_match)  # strip path and extension
      log_by_sample[[sample]] <- character()
    } else if (!is.null(sample)) { # Group the logs under current filename
      log_by_sample[[sample]] <- c(log_by_sample[[sample]], line)
    }
  }

  # log_by_sample is a list of vectors of lines
  too_smalls_by_sample <- sapply(log_by_sample,
                                 function(lines) {
                                   sum(base::grepl("\\[WARN\\].*too small", lines))
                                   })

  # if stranded merge neg and pos (one is actually always 0 but to be sure)
  if (strand != 0) {
    names(too_smalls_by_sample) <- gsub("(_pos|_neg)", "", names(too_smalls_by_sample))
    too_smalls_by_sample <- tapply(too_smalls_by_sample,
                                   names(too_smalls_by_sample),
                                   sum)
  }
  return(too_smalls_by_sample)
}

#' Decompress bedgraph.gz files
#'
#' @param input_dir Path to the folder containing .bedgraph.gz files.
#' @param output_dir Path to the folder where the decompressed .bedgraph files will be saved.
#'                   If NULL (default), decompresses in input_dir.
#' @param grep_pattern Optional. Filter for a grep pattern in files. perl is TRUE.
#'
#' @return Invisibly, a vector of decompressed file paths.
#' @export
decompress_bedgraph_gz <- function(input_dir, output_dir = NULL, grep_pattern = NULL) {
  if (!dir.exists(input_dir)) {
    stop("input_dir does not exist: ", input_dir)
  }
  
  # Default: decompress in place
  if (is.null(output_dir)) {
    output_dir <- input_dir
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  all_files <- list.files(input_dir)
  if (!is.null(grep_pattern)) {
    all_files <- grep(grep_pattern, all_files, perl = TRUE, value = TRUE)
  }
  
  # Find .bedgraph.gz files
  gz_files <- grep("\\.bedgraph\\.gz$", all_files, value = TRUE, ignore.case = TRUE)
  
  if (length(gz_files) == 0) {
    message("No .bedgraph.gz files found.")
    return(invisible(character()))
  }
  
  decompressed_paths <- character()
  
  for (gz_file in gz_files) {
    uncompressed_file <- sub("\\.gz$", "", gz_file)
    full_gz_path <- file.path(input_dir, gz_file)
    full_uncompressed_path <- file.path(output_dir, uncompressed_file)
    
    if (!file.exists(full_uncompressed_path)) {
      message("Decompressing ", gz_file, " to ", output_dir)
      system(paste("gunzip -c", shQuote(full_gz_path), ">", shQuote(full_uncompressed_path)))
    } else {
      message("Skipping ", uncompressed_file, " (already exists)")
    }
    
    decompressed_paths <- c(decompressed_paths, full_uncompressed_path)
  }
  
  invisible(decompressed_paths)
}


# === HELPERS ===

get_first_column <- function(file) {
  #if (endsWith(file, ".gz")) {
  #  cmd <- paste0("zcat ", shQuote(file), " | cut -f1 | uniq")
  #} else {
  #  cmd <- paste0("cut -f1 ", shQuote(file), " | uniq")
  #}
  #return(system(cmd, intern = TRUE))
  unique(data.table::fread(file, select = 1, colClasses = c("character"), 
                           showProgress = TRUE)[[1]])
}

seqname_count_overview <- function(validation_data, anno_seqnames) {
  file_labels <- basename(sapply(validation_data, `[[`, "file"))
  header <- paste(c("seqname", file_labels), collapse = "\t")
  rows <- vapply(anno_seqnames, function(chr) {
    counts <- vapply(validation_data, function(res) {
      v <- res$counts[chr]
      if (is.na(v)) 0L else as.integer(v)
    }, integer(1))
    paste(c(chr, counts), collapse = "\t")
  }, character(1))
  paste(c(header, rows), collapse = "\n")
}

validate_seqnames <- function(annotation, bedgraph_files, logger) {
  # cache for previously validated input combinations
  cache_dir <- tools::R_user_dir("metageneR", which = "cache")

  # Create the directory if it doesn't exist yet
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  input_signature <- list(
    annot_path = annotation,
    annot_mtime = file.info(annotation)$mtime,
    bg_paths   = sort(bedgraph_files),
    bg_mtimes  = file.info(bedgraph_files)$mtime
  )
  key <- digest::digest(input_signature)
  cache_file <- file.path(cache_dir, paste0(key, ".rds"))

  # reprinting messages without reading cols again
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)

    # Invalidate old cache format that lacks counts
    if (!is.list(cached) || is.null(cached$version) || cached$version < 2L) {
      file.remove(cache_file)
    } else {
      logger("[validation seqnames] Skipping check, using cached results.")
      validation_data <- cached$results

      for (res in validation_data) {
        if (!res$is_valid) {
          logger(
            "[validation seqnames] [cached] Not all sequence names in ", res$file, " are present in the annotation.\n",
            "Missing: ", paste(res$missing, collapse = ", "), "\n",
            "Continuing with annotation names...",
            type = "warning"
          )
        }
      }

      logger("[validation seqnames] Row counts per chromosome:\n",
             seqname_count_overview(validation_data, cached$anno_seqnames))
      return(invisible(TRUE))
    }
  }

  # --- VALIDATION LOGIC: COMPUTE AND SAVE ---
  anno_first_col <- get_first_column(annotation)
  results_to_cache <- list()

  for (bedgraph in bedgraph_files) {
    logger("[validation seqnames] Reading the first column of ", bedgraph)
    bg_col1 <- data.table::fread(bedgraph, select = 1, colClasses = c("character"),
                                 showProgress = FALSE)[[1]]
    bg_first_col <- unique(bg_col1)
    bg_counts    <- table(bg_col1)

    missing <- setdiff(bg_first_col, anno_first_col)
    is_valid <- length(missing) == 0

    # Store the state for this specific file
    results_to_cache[[bedgraph]] <- list(
      file     = bedgraph,
      is_valid = is_valid,
      missing  = missing,
      counts   = bg_counts
    )

    if (!is_valid) {
      logger(
        "Not all sequence names in the BEDgraph files are present in the annotation.\n",
        "Missing: ", paste(missing, collapse = ", "), "\n",
        "Continuing with annotation names...",
        type = "warning"
      )
    }
  }

  logger("[validation seqnames] Row counts per chromosome:\n",
         seqname_count_overview(results_to_cache, anno_first_col))
  saveRDS(list(version = 2L, anno_seqnames = anno_first_col, results = results_to_cache),
          file = cache_file)
}

#' Clear the seqname validation cache
#'
#' Deletes all cached seqname validation results written at \code{bin_genome()}
#' step.
#' Use this when annotation or BEDgraph files have been changed.
#'
#' @returns Invisibly returns the number of cache files deleted.
#' @export
clear_validation_cache <- function() {
  cache_dir <- tools::R_user_dir("metageneR", which = "cache")
  if (!dir.exists(cache_dir)) {
    message("No cache directory found: ", cache_dir)
    return(invisible(0L))
  }
  cache_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(cache_files) == 0) {
    message("Cache is already empty: ", cache_dir)
    return(invisible(0L))
  }
  file.remove(cache_files)
  message("Removed ", length(cache_files), " cache file(s) from ", cache_dir)
  invisible(length(cache_files))
}

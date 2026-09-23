.bedtools_env <- new.env(parent = emptyenv())
.bedtools_env$bedtools_path <- "bedtools" # default, used if in PATH

.bw_env <- new.env(parent = emptyenv())
.bw_env$bgtobw_path <- "bigWigToBedGraph" # default, used if in PATH

#' bedtools: utilities for genomic regions
#'
#' @name bedtools
#' @rdname bedtools
NULL

#' @rdname bedtools
#' @section set_bedtools_path:
#'
#' @param path Path to the bedtools executable.
#'
#' @return Returns the path invisibly on successful execution.
#' @export
set_bedtools_path <- function(path) {
  if (!file.exists(path)) stop("Invalid bedtools path: ", path)
  .bedtools_env$bedtools_path <- normalizePath(path)
  invisible(path)
}

#' @rdname bedtools
#' @section bt_genomecov:
#' Executes the following command:
#'  bedtools genomecov -ibam input -bg (-strand strand) > output_file.
#' Creates or returns coverage file based on BAM or BED|GFF|VCF.
#'
#'
#' @param input Path to BAM or BED|GFF|VCF file.
#' @param output_file path to output file.
#' @param stranded logical. Should output be split into strand-specific two files?
#' Uses -strand argument of genomecov.
#' Default is FALSE.
#' @param genome path to genome file. required for non-BAM input, defaults to NULL
#'
#' @details
#' If no bedtools path is given, the executable called is "bedtools" per default and runs if it is in the PATH.
#'
#' The extension of `output_file` will be stripped.
#' If stranded = TRUE, it generates two output files suffixed with "_pos" and "_neg".
#'
#' @returns If output file given, its path or a vector of paths (if stranded) is returned.
#' Otherwise, a data table of the result captured from the command, or
#' if stranded, a list of two data tables: (+,-).
#' @export
#'
#'
#' @examples
#' \dontrun{
#'   set_bedtools_path("/usr/bin/bedtools")
#'   bt_genomecov("bamfile.bam", "output.bedgraph", stranded = TRUE)
#' }
bt_genomecov <- function(input, output_file = NULL, stranded = FALSE, genome = NULL) {
  
  params <- as.list(environment())

  # ================ Config and Logs  ================
  #write_config("bt_genomecov", params)
  #revert_sink <- start_log("bt_genomecov", params) # starts the log and returns a function to revert the sink on exit
  #if (!is.null(revert_sink)) on.exit(revert_sink(), add = T)
  # ================================

  # 0: not stranded, 1: regular, -1: reverse stranded

  # stranded = FALSE: not stranded
  # stranded = TRUE: strand-specific output (split into + and -) -> two files

  bedtools <- .bedtools_env$bedtools_path
  if(grepl("\\.bam$", input, ignore.case = TRUE)) {
    in_param <- c("-ibam", input)
  } else {
    if(is.null(genome)) stop("genome file -g required for non-BAM input")
    in_param <- c("-i", input, "-g", genome)
  }
  base_args <- c("genomecov", in_param, "-bg")

  # Set strand mapping if given
  strand_map <- c(pos = "+", neg = "-")

  if (!is.null(output_file)) {
    output_file <- sub("\\.bedgraph$", "", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

    # if stranded -> split the result into strands
    if (stranded) {

      for (s in names(strand_map)) {
        args <- c(base_args, "-strand", strand_map[[s]])
        file <- paste0(output_file, "_", s, ".bedgraph")
        processx::run(bedtools, args, stdout = file, echo_cmd = T)
      }
    }

    # if not stranded do one call
    else {
      processx::run(bedtools, base_args, stdout = paste0(output_file, ".bedgraph"), echo_cmd = T)
    }

    # return output paths
    if(stranded) {
      return(invisible(paste0(output_file, "_", names(strand_map), ".bedgraph")))
    }
    return(invisible(paste0(output_file, ".bedgraph")))
  }

  # No output file -> return data in R
  if (stranded) {
    result <- lapply(strand_map, function(s) {
      args <- c(base_args, "-strand", s)
      res <- processx::run(bedtools, args)
      data.table::fread(res$stdout)
    })
    names(result) <- names(strand_map)
    return(result)
  }

  res <- processx::run(bedtools, base_args)
  return(data.table::fread(res$stdout))
}

#' @rdname bedtools
#' @section set_bigwigtobedgraph_path:
#'
#' @param path Path to the bigWigToBedGraph executable.
#'
#' @return Returns the path invisibly on successful execution.
#' @export
set_bigwigtobedgraph_path <- function(path) {
  if (!file.exists(path)) stop("Invalid bigWigToBedGraph path: ", path)
  .bw_env$bgtobw_path <- normalizePath(path)
  invisible(path)
}

#' @rdname bedtools
#' @section bw_to_bedgraph:
#' Executes the following command:
#'  bigWigToBedGraph input.bw output.bedgraph
#' Creates or returns bedgraph from BigWig file.
#'
#' @param input Path to BigWig (.bw) file.
#' @param output_file Path to output bedgraph file. If NULL, result returned as data.table.
#' @param chrom Character. Restrict output to given chromosome. Default NULL (all chromosomes).
#' @param start Integer. Start position (0-based) for region restriction. Requires chrom. Default NULL.
#' @param end Integer. End position for region restriction. Requires chrom. Default NULL.
#'
#' @details
#' If no executable path is set, calls "bigWigToBedGraph" and requires it to be in PATH.
#'
#' The extension of `output_file` will be stripped and ".bedgraph" appended.
#'
#' @returns If output_file given, its path is returned invisibly.
#' Otherwise, a data.table of the bedgraph result.
#' @export
#'
#' @examples
#' \dontrun{
#'   set_bigwigtobedgraph_path("/usr/local/bin/bigWigToBedGraph")
#'   bw_to_bedgraph("signal.bw", "output.bedgraph")
#'   bw_to_bedgraph("signal.bw", "output.bedgraph", chrom = "chr1", start = 0, end = 1000000)
#' }
bw_to_bedgraph <- function(input, output_file = NULL, chrom = NULL, start = NULL, end = NULL) {
  stopifnot(endsWith(input, ".bw") || endsWith(input, ".bigwig") || endsWith(input, ".BigWig"))

  if (!is.null(start) || !is.null(end)) {
    if (is.null(chrom)) stop("chrom required when start or end given")
  }

  exe <- .bw_env$bgtobw_path

  build_args <- function(out) {
    args <- input
    if (!is.null(chrom)) args <- c(args, "-chrom", chrom)
    if (!is.null(start)) args <- c(args, "-start", as.character(start))
    if (!is.null(end))   args <- c(args, "-end",   as.character(end))
    c(args, out)
  }

  if (!is.null(output_file)) {
    output_file <- sub("\\.bedgraph$", "", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    out <- paste0(output_file, ".bedgraph")
    processx::run(exe, build_args(out), echo_cmd = TRUE)
    return(invisible(output_file))
  }

  tmp <- tempfile(fileext = ".bedgraph")
  on.exit(unlink(tmp), add = TRUE)
  processx::run(exe, build_args(tmp), echo_cmd = TRUE)
  data.table::fread(tmp)
}

#' @rdname bedtools
#' @section bed_to_bedgraph:
#' Splits a BED with 6 columns (chr, start, end, name, score, strand) with non-overlapping intervals into one bedgraph file per strand.
#' The score from column 5 becomes the bedgraph value. 
#' Easy workaround using neither bedtools nor a genome file.
#' 
#' @param input_file Path to bed file with 6 columns. 
#' @param output_file Path to output file, "_pos.bedgraph" and "_neg.bedgraph" are appended. 
#'  Default NULL: the result is returned as a list of two data.tables pos and neg.
#' 
#' @details
#' Intervals on the same strand must not overlap!
#' Strand column 6 must have values + and -.
#' Score column must be numeric.
#' For overlapping intervals use \code{\link{bt_genomecov()}}
#' 
#' @returns If \code{output_file} is given, the two output paths are returned invisibly. 
#' Otherwise, a list of two data.tables pos and neg is returned.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'   bed_to_bedgraph("regions.bed", "out/regions.bedgraph")
#'   bed_to_bedgraph("regions.bed", "out/regions")
#' }
bed_to_bedgraph <- function(input_file, output_file = NULL) {
  dt <- data.table::fread(input_file, header = FALSE, select = c(1:3, 5:6),
                          col.names = c("chr", "start", "end", "score", "strand"))
  if (!all(dt$strand %in% c("+", "-"))) stop("Strand column (6) must be '+' or '-' in ", input_file)
  if (!is.numeric(dt$score)) stop("Score column (5) is not numeric in ", input_file)
  
  # checking for overlaps --> not allowed here, only possible bt genomecov
  data.table::setorderv(dt, c("strand", "chr", "start"))
  same_group <- dt$chr == data.table::shift(dt$chr) & dt$strand == data.table::shift(dt$strand)
  n_overlap <- sum(same_group & dt$start < data.table::shift(dt$end), na.rm = TRUE)
  if (n_overlap > 0) stop(n_overlap, " overlapping interval(s) in ", input_file, "; use bt_genomecov() instead")
  
  strand_map <- c(pos = "+", neg = "-")
  result <- lapply(strand_map, function(s) dt[dt$strand == s, c("chr", "start", "end", "score"), with = FALSE])
  if (is.null(output_file)) return(result)
  
  output_file <- sub("\\.bedgraph$", "", output_file)
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  out <- paste0(output_file, "_", names(strand_map), ".bedgraph")
  for (i in seq_along(out)) data.table::fwrite(result[[i]], out[i], sep = "\t", col.names = FALSE)
  invisible(out)
}

#' @rdname bedtools
#' @section abs_bedgraph:
#' Converts negative score values in a BEDgraph file to their absolute values.
#'
#' @param input Path to a `.bedgraph` file.
#' @param output_file Path for the output `.bedgraph` file.
#'   If NULL, returns the result as a data.table.
#'
#' @details
#' Only the fourth (score) column is modified. When `output_file` is given, the
#' conversion is done with two streaming `awk` passes — no full file load into R
#' memory. If no negative values are present the file is copied unchanged.
#'
#' @returns If `output_file` is given, its path invisibly. Otherwise a data.table.
#' @export
#'
#' @examples
#' \dontrun{
#'   abs_bedgraph("signal_neg.bedgraph", "signal_abs.bedgraph")
#' }
abs_bedgraph <- function(input, output_file = NULL) {
  stopifnot(endsWith(input, ".bedgraph"))

  if (!is.null(output_file)) {
    output_file <- sub("\\.bedgraph$", "", output_file)
    out <- paste0(output_file, ".bedgraph")
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    count_res <- processx::run("awk", c("$4 < 0 {c++} END{print c+0}", input))
    neg_count <- as.integer(trimws(count_res$stdout))
    if (neg_count == 0) {
      message("[abs_bedgraph] No negative values in ", basename(input))
      file.copy(input, out, overwrite = TRUE)
    } else {
      message("[abs_bedgraph] Converting ", neg_count, " negative value(s) to positive in ", basename(input))
      processx::run("awk", c('BEGIN{OFS="\\t"} {if ($4 < 0) $4 = -$4; print}', input), stdout = out)
    }
    return(invisible(out))
  }

  dt <- data.table::fread(input, header = FALSE)
  neg_count <- sum(dt[[4]] < 0, na.rm = TRUE)
  if (neg_count == 0) {
    message("[abs_bedgraph] No negative values in ", basename(input))
  } else {
    message("[abs_bedgraph] Converting ", neg_count, " negative value(s) to positive in ", basename(input))
    dt[, V4 := abs(V4)]
  }
  dt
}

#' @rdname bedtools
#' @section strip_bedgraph_header:
#' Detects and removes leading comment or header lines from a BEDgraph file.
#'
#' @param input Path to a `.bedgraph` file.
#' @param output_file Path for the output `.bedgraph` file.
#'   If NULL, returns the cleaned result as a data.table.
#'
#' @details
#' Lines at the start of the file beginning with `#`, `track`, or `browser`
#' are considered header/comment lines and are skipped. Their content is printed
#' as an informational message. A file connection reads only the leading lines and
#' stops at the first data line. When `output_file` is given, `tail -n +k` streams
#' the remainder directly to disk without loading the file into R memory.
#'
#' @returns If `output_file` is given, its path invisibly. Otherwise a data.table.
#' @export
#'
#' @examples
#' \dontrun{
#'   strip_bedgraph_header("signal_with_header.bedgraph", "signal_clean.bedgraph")
#' }
strip_bedgraph_header <- function(input, output_file = NULL) {
  stopifnot(endsWith(input, ".bedgraph"))
  con <- file(input, "r")
  on.exit(close(con), add = TRUE)
  header_lines <- character()
  found_data <- FALSE
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (grepl("^\\s*(#|track|browser)", line)) {
      header_lines <- c(header_lines, line)
    } else {
      found_data <- TRUE
      break
    }
  }
  if (!found_data) stop("No data lines found in ", input)
  n_skipped <- length(header_lines)
  if (n_skipped > 0) {
    message("[strip_bedgraph_header] Skipping ", n_skipped, " header line(s) in ", basename(input), ":\n",
            paste(" ", header_lines, collapse = "\n"))
  } else {
    message("[strip_bedgraph_header] No header lines found in ", basename(input))
  }

  if (!is.null(output_file)) {
    output_file <- sub("\\.bedgraph$", "", output_file)
    out <- paste0(output_file, ".bedgraph")
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    processx::run("tail", c(paste0("-n +", n_skipped + 1), input), stdout = out)
    return(invisible(out))
  }

  data.table::fread(input, skip = n_skipped, header = FALSE)
}

# HELPERS =========

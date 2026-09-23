#' Metagene plot
#'
#' @rdname plotting
#' @description
#' Basic plot drawing function wrapped from binGenome plotting library
#'
#' @param plot_dir A character string. Where to save PDF output.
#' @param file_name A character string. File identifier.
#' Currently "metagene_" is automatically added as prefix.
#' @param binMatrixList List of binMatrix instances.
#' @param pvThresholdsVec A vector of numeric between 0 and 1 to use
#' as p-value thresholds in defaultPvalueColorTransformer for Wilcoxon test.
#' @keywords internal
plot_metagene_from_binmatrix <- function(plot_dir, file_name, binMatrixList, pvThresholdsVec = c(0.05, 10^-3, 10^-5), ...) {
  # Check if binMatrixList is named -> necessary since used in the legend
  if (is.null(names(binMatrixList))) stop("binMatrixList should be a named list!")

  # Default arguments
  defaults <- list(
    labels = NULL,
    minYA = 0
  )

  # Unite formal and dot arguments into list
  dots <- list(...)

  pdf_width  <- if (!is.null(dots$width))  dots$width  else 14
  pdf_height <- if (!is.null(dots$height)) dots$height else 7
  pdf_pointsize <- if(!is.null(dots$pointsize)) dots$pointsize else 12
  dots$width  <- NULL
  dots$height <- NULL
  dots$pointsize <- NULL

  # Keep the relevant arguments for wrapped plot function and warn for unused
  all_arguments <- names(formals(plotGroup))
  extra_args <- setdiff(names(dots), all_arguments)
  if (length(extra_args) > 0) {
    warning("Ignoring unused arguments passed to plot_metagene: ", paste(extra_args, collapse = ", "))
  }
  dots <- dots[names(dots) %in% all_arguments]

  # Update arguments by the user's valid arguments
  defaults <- modifyArgs(defaults, dots)

  # Final arguments the call wrapped function
  final_args <- c(mget("binMatrixList", envir = environment()), defaults)

  pdf(paste0(plot_dir, "/metagene_", file_name, ".pdf"), width = pdf_width, height = pdf_height, pointsize = pdf_pointsize)
  on.exit(dev.off(), add = TRUE)

  # plotting.R:51, replacing the if/else
  n_is_two <- length(binMatrixList) == 2
  if (is.null(final_args$performWilcoxTest)) final_args$performWilcoxTest <- n_is_two
  if (final_args$performWilcoxTest && !n_is_two) stop("performWilcoxTest = TRUE needs exactly two conditions.")
  if (final_args$performWilcoxTest) {
    final_args$wilcoxTestPVTransformFUN <- function(x) defaultPvalueColorTransformer(x, pvThresholds = pvThresholdsVec)
  }

  args_to_write <- final_args[names(final_args) != "binMatrixList"]
  args_to_write$filename <- file_name
  args_to_write$binMatrixList_names <- names(binMatrixList)
  jsonlite::write_json(args_to_write,
                       file.path(plot_dir, paste0(file_name, ".json")),
                       pretty = TRUE, auto_unbox = TRUE)

  print_args <- args_to_write
  print_args$binMatrixList <- ""
  do.call(plotGroup, final_args)
}

#' Show parameter options for metagene plotting
#' @rdname plotting
#'
#' @returns Returns no value. Prints the parameters.
#' @export
args_for_metagene_plot <- function() {
  db <- tools::Rd_db("metageneR")
  rd <- db[["plotGroup.Rd"]]
  if (is.null(rd)) return(invisible(NULL))

  # find the \arguments section
  arg_sections <- rd[sapply(rd, function(x) attr(x, "Rd_tag") == "\\arguments")]
  if (length(arg_sections) == 0) return(invisible(NULL))

  items <- arg_sections[[1]]

  # keep only \item elements with a non-empty name
  items <- Filter(function(arg) attr(arg, "Rd_tag") == "\\item" && length(arg[[1]]) > 0, items)

  # skip the first three parameters
  if (length(items) <= 3) return(invisible(NULL))
  items <- items[4:length(items)]

  extract_rd_text <- function(node) {
    if (is.character(node)) return(node)
    if (identical(attr(node, "Rd_tag"), "\\preformatted")) {
      content <- paste(sapply(node, extract_rd_text), collapse = "")
      return(sub("^\n+", "", content))
    }
    paste(sapply(node, extract_rd_text), collapse = "")
  }

  # print remaining args
  for (arg in items) {
    name <- paste(unlist(arg[[1]]), collapse = "")
    desc <- paste(sapply(arg[-1], extract_rd_text), collapse = "")
    desc <- trimws(desc)
    cat(name, ":", desc, "\n")
  }
}

#' Pairwise and all
#'
#' Plots all pairwise comparisons and one with all conditions.
#' Exploits multithreading from package parallel using mclapply on pairs generated.
#' binMatrix objects in the coverage folder in run_dir are aggregated
#' using aggregate_FUN (default: mean) and saved as rds objects in a directory
#' called 'aggregated_binmatrices' as a cache for further calls of the same function
#' with same run_dir, annotation and aggregate function.
#'
#' @param plot_dir Where to save PDF outputs of plotting
#' #' @param plot_prefix Is simply a user-given label for the annotation used for binning.
#' It is not a file path. Appends to the file name (<plot_prefix>_<comparisons>.pdf for single comparisons, all_<plot_prefix>.pdf for all samples in one plot).
#' If a named vector, the name is used for the title!
#' @param run_dir The run folder where binning took place.
#' The folder should have a 'config.json' file where the parameters for bin_genome live.
#' Default place for coverage files is 'out/' but if it does not exist, the function will look for coverage files in the run_dir itself.
#' @param title A character string for the title of the plot.
#' @param metagene_data Output of \code{\link{load_metagene_data}}. Alternative to \code{run_dir}.
#' @param condition_mapper A data frame with columns "grep_name", "display_name", "color" (optional) and "linetype" (optional)
#' Used to match conditions with line colors, line types, and names to display in the legend.
#' If `NULL` (default), conditions are auto-detected either from an existing RDS or coverage (.sum.csv) files of the run folder,
#' and default colors are assigned automatically.
#' This can also be used to specify a subset of conditions to plot. See details.
#' @param aggr_rep_FUN A closure. Function to aggregate the replicates of a sample.
#' @param region_start_label,region_end_label A character string for labels of the region start and end.
#' Default is "TSS" and "TTS", respectively.
#' `region_start_label` is also used as the label for the bin if binning was done with a single bin for the gene body.
#' @param start_ticks,end_ticks An integer. Number of ticks to use for start and end labels if auto-generating (when fixedLabelsStart and fixedLabelsEnd not provided).
#' If `NULL` (default), a suitable default is chosen to evenly divide the upstream and downstream bins, respectively.
#' If the specified value does not evenly divide the bins, the closest integer that evenly divides the bins will be used.
#' @param ids_to_subset A character vector of IDs to filter the regions.
#' @param gsub_name A character string to remove from the filename after parsing. Use for cleaner names in the plot.
#' @param threads Number of threads to use in pairwise plotting.
#' Usage makes sense if multiple pairs are plotted, e.g. nrow(condition_mapper) > 1,
#' and cores are available, e.g. if detectCores() >= threads or
#' upstream threading has not used them up. Since this utilizes mclapply from
#' package parallel, this option should be used careful with nested parallelization.
#' @param overwrite_rds_objects Logical. Should the cached RDS object be overwritten?
#' If `FALSE`, the RDS objects are read, if they exist, to collect the binMatrix objects
#' or are created, if they do not.
#' Default is `FALSE`.
#' @param plot_pairs Controls pairwise comparisons. Use `"all"` to plot all pairwise combinations,
#' or a list of 2-element character vectors to specify explicit pairs
#' (e.g. `list(c("WT", "KO"), c("WT", "DKO"))`). Each element must match a `display_name` in
#' the condition_mapper (equivalently, a name in the subsetted RDS object).
#' If not provided (default), no pairwise plots are generated.
#' @param upstream_bp,downstream_bp Numeric. Lengths in base pairs for the upstream and downstream regions to be shown in the plot, respectively.
#' This is useful if the feature's real upstream and downstream regions differ from the ones specified in the binning procedure.
#' If `NULL` (default), the full upstream/downstream length from the binning config is used automatically.
#' Ignored if `fixedLabelsStart` and `fixedLabelsEnd` are provided, respectively.
#' @param cov_files_to_ignore A regex string passed to `grep(..., invert = TRUE)` to exclude matching coverage files from being read.
#' Default is `".norm.coverage.csv"` to skip normalized coverage files when raw and normalized files coexist in the coverage folder.
#' Set to `NULL` to disable filtering.
#'   When provided, \code{condition_mapper}, \code{aggr_rep_FUN}, \code{ids_to_subset},
#'   \code{gsub_name}, \code{threads}, \code{overwrite_rds_objects}, and
#'   \code{cov_files_to_ignore} are ignored.
#' @param normByShapeSum Logical. Normalize by shape sum (intensity). Default \code{FALSE}.
#' @param normByShapeMax Logical. Normalize by shape max. Default \code{FALSE}.
#' @param aggregateFun Function to aggregate across features per bin. Default \code{mean}.
#'  The y-axis label is derived automatically for \code{mean} and \code{median}, for any other function a generic label is used. PAss \code{ylab} to set a custom label.
#' @param ... Further arguments for plotting.
#' Call `args_for_metagene_plot()` to see the options.
#'
#' @details The function relies on the structure of the run directory and the config file to understand where the coverage files are and how to plot.
#' The `run_dir` should have a `config.json` file where the parameters for bin_genome live.
#' The default place for coverage files is `out/` but if it does not exist, the function will look for coverage files in the `run_dir` itself.
#'
#' The list of grouped binMatrix objects (that is cached by saving as RDS) contains the aggregates across replicates.
#' The names of this list correspond to the `display_name` in the `condition_mapper`.
#' Subsetting is done by the display name using an exact match (no grep or regex), to avoid susbtring matches and allow for use of special characters.
#'
#' condition_mapper is automatically generated if not provided by the user.
#' For this purpose, the basename of `.sum.csv` files in the coverage folder are used to extract the unique condition names.
#'
#' Both at the annotation step and the binning step, upstream and downstream regions are defined:
#' `upstream,downstream` parameters in `make_windows()` and fixedBinSize parameters in `bin_genome()`.
#' This allows the start and end label regions in the plot to extend into the gene body.
#' For this purpose, binning is done by longer fixed bin sizes than the "real" upstream and downstream regions specified for the annotation.
#' Alternatively, the plotting parameters `fixedLabelsStart` and `fixedLabelsEnd` specify
#' the labels for upstream and downstream (see `args_for_metagene_plot()` for more details).
#' The last element of fixedLabelsStart and the first element of fixedLabelsEnd overwrite the region_start_label and region_end_label, respectively.
#'
#' @returns No return value. Creates plots in the `plot_dir`.
#'
#' @export
plot_metagene_experiment <- function(plot_dir, plot_prefix, run_dir = NULL, title = plot_prefix,
                                     metagene_data = NULL, 
                                     condition_mapper = NULL, aggr_rep_FUN = mean,
                                     region_start_label = "TSS", region_end_label = "TTS",
                                     start_ticks = NULL, end_ticks = NULL,
                                     ids_to_subset = NULL, gsub_name = NULL,
                                     threads = 1,
                                     overwrite_rds_objects = FALSE,
                                     plot_pairs = NULL, upstream_bp = NULL, downstream_bp = NULL,
                                     cov_files_to_ignore = ".norm.coverage.csv",
                                     normByShapeSum = FALSE, normByShapeMax = FALSE,
                                     aggregateFun = mean,
                                     ...) {
  
  # check arguments
  if(is.null(run_dir) && is.null(metagene_data)) stop("Provide either run_dir or metagene_data.")
  if(!is.null(run_dir) && !is.null(metagene_data)) stop("Provide either run_dir or metagene_data, not both.")
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  
  # load data
  if(is.null(metagene_data)) {
    metagene_data <- load_metagene_data(
      run_dir = run_dir, condition_mapper = condition_mapper, aggr_rep_FUN = aggr_rep_FUN,
      ids_to_subset = ids_to_subset, gsub_name = gsub_name, threads = threads,
      overwrite_rds_objects = overwrite_rds_objects, cov_files_to_ignore = cov_files_to_ignore
    )
  }
  grouped <- metagene_data
  config <- attr(metagene_data, "config")
  condition_mapper <- attr(metagene_data, "condition_mapper")
  condition_greps <- condition_mapper[["grep_name"]]
  
  if(is.null(config)) stop("metagene_data has no config attribute. Was it created by load_metagene_data()?")
  
  # validate plot_pairs
  if(!is.null(plot_pairs)) {
    if(identical(plot_pairs, "all")) {
      message("[plot_metagene_experiment] comparing all pairs")
    } else if(is.list(plot_pairs)) {
      bad <- which(!vapply(plot_pairs, function(p) is.character(p) && length(p) == 2L, logical(1)))
      if(length(bad) > 0) stop("[plot_metagene_experiment] plot_pairs elements at index", paste(bad, collapse = ", "))
    } else {
      stop("[plot_metagene_experiment] plot_pairs must be NULL, all, or a list of 2-element character vectors.")
    }
  }
  
  if(is.list(plot_pairs)) {
    valid_names <- condition_mapper$display_name
    for(i in seq_along(plot_pairs)) {
      unknown <- setdiff(plot_pairs[[i]], valid_names)
      if(length(unknown) > 0) stop("[plot_metagene_experiment] plot_pairs[[", i, "]] contains unknown condition(s): ", 
                                   paste(unknown, collapse = ", "), ".\n Please use ", paste(valid_names, collapse = ", "))
    }
  }
  
  dotArgs <- list(...)
  if (any(c("binMatrixList", "labels", "title", "normByLibSize", "normByBinlength") %in% names(dotArgs))) {
    warning("[plot_metagene_experiment] Arguments binMatrixList, labels, and title as well as normByLibSize and normByBinlength cannot be set via ... Ignoring.")
    dotArgs <- dotArgs[!names(dotArgs) %in% c("binMatrixList", "labels", "title", "normByLibSize", "normByBinlength")]
  }
  
  message("[plot_metagene_experiment] Started processing with annotation: ", plot_prefix)
  message("[plot_metagene_experiment] Using the following conditions: ", paste(condition_mapper$display_name, collapse = ", "))

  # === label settings ===
  if (any(c("fixedLabelsStartTotalBins", "fixedLabelsEndTotalBins") %in% names(dotArgs))) {
    warning("[plot_metagene_experiment] You provided `fixedLabelsStartTotalBins` or `fixedLabelsEndTotalBins` which is going to be ignored, since it is only read from the config.")
  }
  axis_labels <- metagene_axis_labels(config, region_start_label, region_end_label, start_ticks, end_ticks,
                                      upstream_bp, downstream_bp, dotArgs[["fixedLabelsStart"]], dotArgs[["fixedLabelsEnd"]])

  # === Plotting ===
  message("[plot_metagene_experiment] Plotting with condition_mapper: ")
  print(condition_mapper)

  # Update dotArgs with computed local variables so they override user-provided values
  dotArgs <- modifyArgs(dotArgs, c(axis_labels, list(
    title = title,
    normByShapeSum = normByShapeSum,
    normByShapeMax = normByShapeMax,
    aggregateFun = aggregateFun
  )))

  # Pairwise plotting if wanted
  if (!is.null(plot_pairs)) {
    if (identical(plot_pairs, "all")) {
      # sanity check that there are actually more than two conds
      if(length(condition_greps) < 2) {
        warning("[plot_metagene_experiment] plot_pairs = 'all' needs at least 2 conditions, skipping pairwise plots")
      } else {
        pairs_matrix <- combn(condition_greps, 2)
        pairs_list <- lapply(seq_len(ncol(pairs_matrix)), function(i) {
          greps <- pairs_matrix[, i]
          idx <- match(greps, condition_mapper$grep_name)
          condition_mapper$display_name[idx]
        })
      }
    } else {
      pairs_list <- plot_pairs
    }

    message("[plot_metagene_experiment] Pairwise plotting in threads: ", threads)
    parallel_applier(seq_along(pairs_list), function(comparison_index) {
      display_names <- pairs_list[[comparison_index]]

      color_idx <- match(display_names, condition_mapper$display_name)
      colors <- condition_mapper$color[color_idx]
      linetypes <- condition_mapper$linetype[color_idx]

      message("[plot_metagene_experiment] Drawing plot for ", paste(display_names, collapse = ", "),
              " with color(s) ", paste(colors, collapse = ", "), ".")

      pairwise_comparison <- grouped[display_names]

      do.call(plot_metagene_from_binmatrix, c(list(
        plot_dir = plot_dir,
        file_name = paste0(plot_prefix, "_", paste0(gsub(" ", "_", display_names), collapse = "_vs_")),
        binMatrixList = pairwise_comparison,
        palette = colors,
        shapeLineType = linetypes
      ), dotArgs))

      return(paste(display_names, collapse = "_vs_"))
    }, cores = threads, preschedule = TRUE)
  }

  # Plot all genes:
  do.call(plot_metagene_from_binmatrix, c(list(
    plot_dir = plot_dir,
    file_name = paste0("all_", plot_prefix),
    binMatrixList = grouped,
    palette = condition_mapper$color,
    shapeLineType = condition_mapper$linetype
  ), dotArgs))
}


#' Re-plot using pars in json
#'
#' This function enables re-plotting of a metagene plot
#' by using parameters given in a JSON file and by specifying parameters
#' to overwrite those in the JSON.
#' The JSON file should have an attribute called binMatrixList_names, which
#' indicate the conditions to be compared.
#'
#' @param plot_dir Where to save PDF of the plot
#' @param file_name File identifier. "metagene_" gets added as prefix.
#' @param path_to_binMatrix Path to the RDS object of aggregated binMatrices
#' that has been created by `plot_metagene_experiment()`.
#' @param json_file JSON file with a list of parameters to use in plotting
#' @param ... Parameters to set/overwrite those in the JSON. Same as `plot_metagene_experiment()`.
#'
#' @export
plot_metagene_from_json <- function(plot_dir, file_name, path_to_binMatrixList, json_file, ids_to_subset = NULL, ...) {
  # Argument validity
  if (!file.exists(path_to_binMatrixList)) stop("Path to the binMatrix RDS object not found: ", path_to_binMatrixList)
  if (!grepl("\\.(R|r)ds$", path_to_binMatrixList)) stop("Not an RDS file: ", path_to_binMatrixList)
  if (!file.exists(json_file)) stop("JSON file not found: ", json_file)

  # Read RDS object created by plot_metagene_experiment
  message("Reading RDS object: ", path_to_binMatrixList)
  grouped <- readRDS(path_to_binMatrixList)

  if (!is.null(ids_to_subset)) {
    # filterByIDs works on binMatrix
    # since grouped has conditions with *aggregated* bin matrices,
    # iterate over conditions
    grouped <- lapply(grouped, filterByIDs, filterIDs = ids_to_subset)
    gc() # lower memory overhead
  }
  
  # subsetting stats (how many rows from binmatrix are saved, how many ids from subsetting vector are actually used)
  if(!is.null(ids_to_subset)) {
    nr_requested <- length(unique(ids_to_subset))
    ids_kept <- rownames(grouped[[1]]@cov)
    nr_missing <- length(setdiff(unique(ids_to_subset), ids_kept))
    if(length(ids_kept) == 0) stop("[plot_metagene_from_json] None of the requested ids in <ids_to_subset> are in the binMatrices, check format.")
    message("[plot_metagene_from_json] ", nr_requested, " ids_to_subset -> ", length(ids_kept), " ids kept in binMatrices + ", nr_missing, " ids not present in the data")
  }

  # Parse parameters from JSON
  plot_parameters <- jsonlite::read_json(json_file, simplifyVector = TRUE)
  # Fetch filename from the json if not given and exists in json
  if (is.null(file_name) && "filename" %in% names(plot_parameters)) file_name <- plot_parameters$filename

  # Check if the names of RDS match with names specified in the JSON
  if (!("binMatrixList_names" %in% names(plot_parameters))) stop("JSON file does not have binMatrixList_names attribute!")
  binMatrixList_names <- plot_parameters$binMatrixList_names
  if (!all(binMatrixList_names %in% names(grouped))) stop("binMatrixList_names of JSON file does not comply with RDS binMatrixList. The RDS object has names: ", paste(names(grouped)))
  pairwise_comparison <- grouped[binMatrixList_names]

  # remove names attribute unnecessary for plotting
  plot_parameters$binMatrixList_names <- NULL

  # === Overwrite Parameters with Dots === #
  dots <- list(...)
  plot_parameters <- modifyArgs(plot_parameters, dots)

  all_args <- c(list(plot_dir = plot_dir, file_name = file_name, binMatrixList = pairwise_comparison),
                plot_parameters)

  do.call(plot_metagene_from_binmatrix, all_args)
}

#' metagene plots per feature cluster
#' 
#' generates a pdf with one metagene panel per id set. Metagene panels contain all given conditions and have the same y axis. 
#' Panel titles are the names of \code{id_sets} and the nr of features in each cluster. 
#' Other behaviour is like \code{\link{plot_metagene_experiment}}, refer to this documentation for more information.
#' 
#' @inheritParams plot_metagene_experiment
#' @param file_name Character string, appended to \code{"metagene_clusters_"} for the PDF name.
#' @param id_sets Named list of character verctors conaining feature ids, possible source \code{$clusters} from \code{\link{cluster_metagene}}
#' @param pvThresholdsVec Numeric vetor of p-value thresholds for colors of the wilcoxon strip (for 2 conditions only)
#' @param ... Further plotting arguments, see \code{args_for_metagene_plot()}. \code{width} and \code{height} set the PDF size (default 14 x 4 per panel, pass
#'  \code{height} for many sets). Panel titles are the names of \code{id_sets}; \code{title} is not supported. The legend is drawn in the
#'  first panel and the x-axis in the last one; pass \code{showLegend}/\code{showXLab} (one value per panel)
#'  to change that.
#' @export
plot_metagene_clusters <- function(plot_dir, file_name, id_sets, metagene_data = NULL, run_dir = NULL, 
                                   condition_mapper = NULL, aggr_rep_FUN = mean, region_start_label = "TSS", region_end_label = "TTS", 
                                   start_ticks = NULL, end_ticks = NULL, upstream_bp = NULL, downstream_bp = NULL, 
                                   ids_to_subset = NULL, gsub_name = NULL, threads = 1, 
                                   overwrite_rds_objects = FALSE, cov_files_to_ignore = ".norm.coverage.csv", 
                                   normByShapeSum = FALSE, normByShapeMax = FALSE, aggregateFun = mean, 
                                   pvThresholdsVec = c(0.05, 10^-3, 10^-5), ...) {
  
  # check arguments
  if(is.null(run_dir) && is.null(metagene_data)) stop("Provide either run_dir or metagene_data.")
  if(!is.null(run_dir) && !is.null(metagene_data)) stop("Provide either run_dir or metagene_data, not both.")
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  
  # load data
  if(is.null(metagene_data)) {
    metagene_data <- load_metagene_data(
      run_dir = run_dir, condition_mapper = condition_mapper, aggr_rep_FUN = aggr_rep_FUN,
      ids_to_subset = ids_to_subset, gsub_name = gsub_name, threads = threads,
      overwrite_rds_objects = overwrite_rds_objects, cov_files_to_ignore = cov_files_to_ignore
    )
  }
  config <- attr(metagene_data, "config")
  condition_mapper <- attr(metagene_data, "condition_mapper")
  if(is.null(config)) stop("metagene_data has no config attribute. Was it created by load_metagene_data()?")
  
  ids_present <- rownames(metagene_data[[1]]@cov)
  for(set in names(id_sets)) {
    nr_found <- sum(unique(id_sets[[set]]) %in% ids_present)
    message("[plot_metagene_clusters] ", set, " has ", nr_found, " of ", length(unique(id_sets[[set]])), " ids in the data")
    if(nr_found == 0) stop("[plot_metagene_clusters] none of the ids of set ", set, " are present in the matrix.")
  }
  
  dotArgs <- list(...)
  pdf_width  <- if(!is.null(dotArgs$width))  dotArgs$width  else 14
  pdf_height <- if(!is.null(dotArgs$height)) dotArgs$height else 4 * length(id_sets)
  dotArgs$width  <- NULL
  dotArgs$height <- NULL
  reserved <- c("binMatrixList", "interestLists", "bodyLabels", "fixedLabelsStartTotalBins", "fixedLabelsEndTotalBins",
                "palette", "shapeLineType", "performWilcoxTest", "wilcoxTestPVTransformFUN",
                "normByLibSize", "normByBinlength", "name", "asIsYlab",
                "metaDataNorm", "isListGeneList", "transGeneMapping", "maxTranscripts", "removeEndings", "psCount", "normName", "applyAggregatedNorm")
  bad_params <- intersect(names(dotArgs), reserved)
  if(length(bad_params) > 0) {
    warning("[plot_metagene_clusters] Ignoring args that are not supported or set from config: ", paste(bad_params, collapse = ","))
    dotArgs <- dotArgs[!names(dotArgs) %in% bad_params]
  }
  
  axis_labels <- metagene_axis_labels(config, region_start_label, region_end_label, start_ticks, end_ticks,
                                      upstream_bp, downstream_bp, dotArgs[["fixedLabelsStart"]], dotArgs[["fixedLabelsEnd"]])
  
  n_panels <- length(id_sets)
  args <- list(showLegend = c(TRUE, rep(FALSE, n_panels - 1)),
               showXLab = c(rep(FALSE, n_panels - 1), TRUE),
               showYLab = TRUE)
  args <- modifyArgs(args, dotArgs)
  args <- modifyArgs(args, list(
    binMatrixList = metagene_data,
    interestLists = id_sets,
    bodyLabels = axis_labels$labels,
    fixedLabelsStart = axis_labels$fixedLabelsStart,
    fixedLabelsStartTotalBins = axis_labels$fixedLabelsStartTotalBins,
    fixedLabelsEnd = axis_labels$fixedLabelsEnd,
    fixedLabelsEndTotalBins = axis_labels$fixedLabelsEndTotalBins,
    palette = condition_mapper$color,
    shapeLineType = condition_mapper$linetype,
    normByShapeSum = normByShapeSum,
    normByShapeMax = normByShapeMax,
    aggregateFun = aggregateFun,
    performWilcoxTest = length(metagene_data) == 2,
    wilcoxTestPVTransformFUN = function(x) defaultPvalueColorTransformer(x, pvThresholds = pvThresholdsVec)
  ))
  
  pdf(file.path(plot_dir, paste0("metagene_clusters_", file_name, ".pdf")), width = pdf_width, height = pdf_height)
  on.exit(dev.off(), add = TRUE)
  do.call(plotMergedShape, args)
  invisible(NULL) # suppresses ugly print of plotMergedShape
}

#' Load metagene data from a run directory
#'
#' Reads binning results from a run directory, groups replicates by condition,
#' and returns the aggregated data. The result can be passed to
#' \code{\link{plot_metagene_experiment}} or \code{\link{get_metagene_profiles}}
#' to avoid re-loading from disk.
#'
#' @param run_dir The run folder where binning took place.
#'   Must contain a \code{config.json} file.
#' @param condition_mapper A data frame with columns \code{grep_name}, \code{display_name},
#'   and optionally \code{color} and \code{linetype}. If \code{NULL} (default),
#'   conditions are auto-detected from coverage files or a cached RDS.
#' @param aggr_rep_FUN Function to aggregate replicates. Default \code{mean}.
#' @param ids_to_subset Character vector of IDs to filter regions.
#' @param gsub_name Character string to remove from filenames for cleaner names.
#' @param threads Number of threads for parallel processing.
#' @param overwrite_rds_objects Logical. Overwrite cached RDS? Default \code{FALSE}.
#' @param cov_files_to_ignore Regex to exclude coverage files. Default \code{".norm.coverage.csv"}.
#'
#' @returns A named list of aggregated binMatrix objects (one per condition).
#'   Can be passed as \code{metagene_data} to \code{plot_metagene_experiment}
#'   or \code{get_metagene_profiles}.
#'
#' @export
load_metagene_data <- function(run_dir, condition_mapper = NULL, aggr_rep_FUN = mean, ids_to_subset = NULL, 
                               gsub_name = NULL, threads = 1, overwrite_rds_objects = FALSE, 
                               cov_files_to_ignore = ".norm.coverage.csv") {
  # sanity checks
  if(!dir.exists(run_dir)) stop("Run dir does not exist: ", run_dir)
  if(!file.exists(file.path(run_dir, "config.json"))) stop("config.json does not exist in run dir")
  
  config <- jsonlite::read_json(file.path(run_dir, "config.json"), simplifyVector = TRUE)
  
  cov_folder <- NULL
  if (dir.exists(file.path(run_dir, "out"))) {
    cov_folder <- file.path(run_dir, "out")
  } else if (dir.exists(run_dir)) {
    cov_folder <- run_dir
  } else {
    stop("Neither 'out/' nor the run_dir itself exists for coverage files. Check your run_dir: ", run_dir)
  }
  cov_files <- list.files(path = cov_folder, pattern = "\\.sum\\.csv$", full.names = FALSE)
  # check if there are coverage files
  if (length(cov_files) == 0) {
    stop("No coverage files found in the coverage folder: ", cov_folder)
  }
  
  # Cached RDS path -> use to auto-generate condition_mapper and for fast reading
  aggr_func <- getFuncName(aggr_rep_FUN)
  rds_object <- file.path(run_dir, paste0("matrix_aggr_with_", aggr_func, ".Rds"))
  
  if(is.null(condition_mapper)){
    message("[load_metagene_data] Autogenerating condition_mapper.")
    if(file.exists(rds_object)) {
      message("[load_metagene_data] Reading condition names from existing RDS: ", rds_object)
      grouped_temp <- readRDS(rds_object)
      conditions <- names(grouped_temp)
    } else {
      message("[load_metagene_data] Extracting condition names from coverage files in ", cov_folder)
      conditions <- unique(gsub("\\.sum\\.csv$", "", basename(cov_files)))
      conditions <- sort(conditions)
    }
    condition_mapper <- data.frame(
      grep_name = conditions, 
      color = grDevices::hcl.colors(length(conditions), palette = "Set 2"), 
      linetype = rep(1, length(conditions)), 
      display_name = conditions, 
      stringsAsFactors = FALSE)
  }
  
  # Validate condition_mapper
  required_columns <- c("grep_name", "display_name")
  absent_columns_in_mapper <- setdiff(required_columns, names(condition_mapper))
  
  if(length(absent_columns_in_mapper) > 0)
    stop("condition_mapper data frame does not have the required columns: ", paste(absent_columns_in_mapper, collapse = ", "))
  
  if(!("color" %in% names(condition_mapper))) {
    condition_mapper$color <- grDevices::hcl.colors(nrow(condition_mapper), palette = "Set 2")
  }
  
  if(!("linetype" %in% names(condition_mapper))) {
    condition_mapper$linetype <- rep(1, nrow(condition_mapper))
  }
  
  condition_greps <- condition_mapper[["grep_name"]]
  
  if(!file.exists(rds_object) || overwrite_rds_objects) {
    message("[load_metagene_data] No RDS object found or overwrite-mode on: ", rds_object)
    grouped <- create_binmatrix(cov_folder, condition_greps, condition_mapper, 
                                aggregate_FUN = aggr_rep_FUN, ids_to_subset = ids_to_subset, 
                                gsub_name = gsub_name, ignore = cov_files_to_ignore, cores = threads)
    if (is.null(ids_to_subset) && is.null(gsub_name)) {
      # only saving if run dir is writable
      cached <- tryCatch({ 
        suppressWarnings(saveRDS(grouped, file = rds_object)); TRUE 
        }, error = function(e) FALSE)
      if (cached) {
        message("[load_metagene_data] Cached RDS for next time: ", rds_object)
      } else {
        message("[load_metagene_data] Not caching, run_dir is not writable: ", run_dir)
      }
    } else {
      message("[load_metagene_data] Not caching, data was created with ids_to_subset/gsub_name; the RDS cache must hold the full data.")
    }
  } else {
    message("[load_metagene_data] Reading RDS object: ", rds_object)
    grouped <- readRDS(rds_object)
    
    # Subset RDS based on samples given in $display_name of condition_mapper
    missing_from_rds <- setdiff(condition_mapper$display_name, names(grouped))
    if (length(missing_from_rds) > 0) {
      in_rds_not_requested <- setdiff(names(grouped), condition_mapper$display_name)
      stop("[load_metagene_data] Some requested conditions not found in RDS.\n",
           "    requested but not in RDS: ", paste(missing_from_rds, collapse = ", "), "\n",
           "    In RDS but not requested: ", paste(in_rds_not_requested, collapse = ", "))
    }
    to_drop <- setdiff(names(grouped), condition_mapper$display_name)
    grouped[to_drop] <- NULL
    grouped <- grouped[condition_mapper$display_name]

    if (!is.null(ids_to_subset)) {
      # filterByIDs works on binMatrix
      # since grouped has conditions with *aggregated* bin matrices,
      # iterate over conditions
      grouped <- lapply(grouped, filterByIDs, filterIDs = ids_to_subset)
      gc() # lower memory overhead
    }
  }
  
  # subsetting stats (how many rows from binmatrix are saved, how many ids from subsetting vector are actually used)
  if(!is.null(ids_to_subset)) {
    nr_requested <- length(unique(ids_to_subset))
    ids_kept <- rownames(grouped[[1]]@cov)
    nr_missing <- length(setdiff(unique(ids_to_subset), ids_kept))
    if(length(ids_kept) == 0) stop("[load_metagene_data] None of the requested ids in <ids_to_subset> are in the binMatrices, check format.")
    message("[load_metagene_data] ", nr_requested, " ids_to_subset -> ", length(ids_kept), " ids kept in binMatrices + ", nr_missing, " ids not present in the data")
  }
  
  # sanity check after subsetting 
  for(condition in grouped) {
    if(is.null(rownames(condition@cov))) stop("Rownames of aggregated bin table are NULL, did you subset correctly?")
  }
  
  # add config as hidden metadata
  attr(grouped, "config") <- config
  attr(grouped, "condition_mapper") <- condition_mapper
  return(grouped)
}

#' Extract metagene profiles as a data.frame
#'
#' Returns the aggregated per-bin signal profiles that would be plotted
#' by \code{\link{plot_metagene_experiment}}, as a data.frame for downstream
#' analysis or custom plotting (e.g. with ggplot2).
#'
#' @param run_dir The run folder where binning took place. Provide either
#'   \code{run_dir} or \code{metagene_data}, not both.
#' @param metagene_data Output of \code{\link{load_metagene_data}}. Provide either
#'   this or \code{run_dir}.
#' @param condition_mapper See \code{\link{plot_metagene_experiment}}.
#'   Ignored when \code{metagene_data} is provided.
#' @param aggr_rep_FUN Function to aggregate replicates. Default \code{mean}.
#'   Ignored when \code{metagene_data} is provided.
#' @param ids_to_subset Character vector of IDs to filter regions.
#'   Ignored when \code{metagene_data} is provided.
#' @param gsub_name Character string to remove from filenames.
#'   Ignored when \code{metagene_data} is provided.
#' @param threads Number of threads. Ignored when \code{metagene_data} is provided.
#' @param overwrite_rds_objects Logical. Ignored when \code{metagene_data} is provided.
#' @param cov_files_to_ignore Regex to exclude files. Ignored when \code{metagene_data} is provided.
#' @param normByShapeSum Logical. Normalize by shape sum. Default \code{FALSE}.
#' @param normByShapeMax Logical. Normalize by shape max. Default \code{FALSE}.
#' @param aggregateFun Function to aggregate across features per bin. Default \code{mean}.
#' @param format Output format: \code{"wide"} (default, one column per condition)
#'   or \code{"long"} (columns: condition, bin, value).
#'
#' @returns A data.frame. Wide: columns \code{bin}, plus one per condition.
#'   Long: columns \code{condition}, \code{bin}, \code{value}.
#'
#' @export
get_metagene_profiles <- function(run_dir = NULL, metagene_data = NULL, condition_mapper = NULL, 
                                  aggr_rep_FUN = mean, ids_to_subset = NULL, gsub_name = NULL, 
                                  threads = 1, overwrite_rds_objects = FALSE, cov_files_to_ignore = ".norm.coverage.csv", 
                                  normByShapeSum = FALSE, normByShapeMax = FALSE,
                                  aggregateFun = mean, format = c("wide", "long")) {
  format <- match.arg(format)
  if(is.null(run_dir) && is.null(metagene_data)) stop("Provide either run_dir or metagene_data.")
  if(!is.null(run_dir) && !is.null(metagene_data)) stop("Provide either run_dir or metagene_data, not both.")
  
  if(!is.null(run_dir)) {
    metagene_data <- load_metagene_data(run_dir = run_dir, condition_mapper = condition_mapper,
                                        aggr_rep_FUN = aggr_rep_FUN, ids_to_subset = ids_to_subset,
                                        gsub_name = gsub_name, threads = threads,
                                        overwrite_rds_objects = overwrite_rds_objects,
                                        cov_files_to_ignore = cov_files_to_ignore)
  }
  
  profiles <- lapply(names(metagene_data), function(cond) {
    meta <- getSummedShape(metagene_data[[cond]], normByShapeSum = normByShapeSum, normByShapeMax = normByShapeMax, aggregateFun = aggregateFun)
    meta@attributes[["meta"]]
  })
  names(profiles) <- names(metagene_data)
  
  if(format == "wide") {
    out <- as.data.frame(profiles)
    out <- cbind(bin = seq_len(nrow(out)), out)
  } else {
    n_bins <- length(profiles[[1]])
    out <- data.frame(
      condition = rep(names(profiles), each = n_bins), 
      bin = rep(seq_len(n_bins), times = length(profiles)), 
      value = unlist(profiles, use.names = FALSE), 
      stringsAsFactors = FALSE
    )
  }
  return(out)
}

# ===== HELPERS =====

#' Collect binTables
#'
#' @param folder Path to the directory of coverage files.
#' @param condition A character string to filter the condition/sample.
#' @param annotation A character string to filter the annotation/region.
#' @param gsub_name A character string to remove from the filename after parsing.
#' @param ids_to_subset A character vector to subset regions
#' @param ignore A character string to filter out files. Generally used for .norm files.
#'
#' @returns a list of binTables.
getBinTable <- function(folder, condition, gsub_name = NULL, ids_to_subset = NULL, ignore = NULL, cores = 1){
  # Does the folder exist?
  if(!dir.exists(folder)) stop("Coverage folder for binTables does not exist! Check: ", folder)

  # filter with condition in the name
  rep_list <- list.files(path = folder,
                         pattern = paste0(condition, ".*\\.coverage\\.csv"),
                         full.names = T)

  # filter out the files with the name to ignore
  if(!is.null(ignore)) rep_list <- grep(ignore, rep_list, invert = T, value = TRUE)

  # fail if empty
  if (length(rep_list) == 0) {
    stop(sprintf("[getBinTable] No files found for condition '%s'. Check your spelling or folder path.", condition))
  }

  message("[getBinTable] Replicates found for condition '", condition, "': ", paste(basename(rep_list), collapse = ", "))

  # filter sum files as well
  cov_list <- list.files(path = folder,
                         pattern = paste0(condition, ".*\\.sum\\.csv"),
                         full.names = T)

  # Check that we have matching coverage and sum files
  if (length(rep_list) != length(cov_list)) {
    stop("[getBinTable] Mismatch: found ", length(rep_list), " coverage files but ", length(cov_list),
         " sum files for condition '", condition, "'. ",
         "Coverage files: ", paste(basename(rep_list), collapse = ", "), ". ",
         "Sum files: ", paste(basename(cov_list), collapse = ", "))
  }

  # Create names from coverage file basenames
  cond_names <- if (is.null(gsub_name)) basename(rep_list) else gsub(gsub_name, "", basename(rep_list))

  if(length(cond_names) == 0) stop("[getBinTable] Names of binTables are empty! Did you mistype the condition? Condition: ", condition)

  names(rep_list) <- cond_names
  names(cov_list) <- cond_names

  message(sprintf("[getBinTable] Reading %d bin tables using %d cores...", length(rep_list), cores))

  cond_list <- parallel_applier(cond_names, function(x) {
    message("[getBinTable]     Bin table name: ", x,
            "\n                     Reading replicate ", rep_list[[x]],
            "\n                     with the coverage sum: ", cov_list[[x]])

    res <- readBinTable(rep_list[[x]], cov_list[[x]], ids_to_subset = ids_to_subset)

    gc()
    return(res)
    }, cores = cores)

  # Are bin tables actually empty? (should not happen)
  if(length(cond_list) == 0) stop("List of binTables is found to be empty! Did you mistype the condition? Condition: ", condition)

  names(cond_list) <- cond_names
  return(cond_list)
}



#' Create binMatrix objects
#'
#' Read from the coverage folder,
#' group by conditions and
#' aggregate results using the provided aggregation function.
#'
#' @param cov_folder Coverage folder, where the binning results are
#' @param condition_greps Patterns for each condition (will match coverage file names)
#' @param annotation_grep A character string to filter the annotation/region.
#' @param condition_mapper See plot_metagene_experiment
#' @param aggregate_FUN See plot_metagene_experiment
#'
#' @returns a list of aggregated binMatrix table results
create_binmatrix <- function(cov_folder, condition_greps, condition_mapper, aggregate_FUN,
                             ids_to_subset = NULL, gsub_name = NULL, ignore = ".norm.coverage.csv", cores = 1) {
  # read bin tables from the coverage folder
  message("Cores specified for BinMatrix creation: ", cores)
  bin_tables <- parallel_applier(condition_greps, function(cond_name) {
    res <- getBinTable(folder = cov_folder,
                condition = cond_name,
                ignore = ignore, # TODO: hard-coded ignore
                gsub_name = gsub_name,
                ids_to_subset = ids_to_subset,
                cores = cores)
    gc()
    return(res)
  }, cores = 1) # no parallelism, because this might cause out-of-memory

  # Assign names for named access: bin_tables[[display_name]]
  names(bin_tables) <- condition_mapper$display_name

  gc() # garbage collection to get some space in case of memory shortage

  # summarize each condition using aggregate_FUN on replicates
  grouped <- parallel_applier(bin_tables, function(condition) {
    message("=== Aggregating bin tables for ", paste(names(condition), collapse = ", "))
    groupMatrices(condition, aggregate_FUN)
    }, cores = cores)

  gc() # garbage collection to remove bin tables
  return(grouped)
}

parallel_applier <- function(X, FUN, cores = 1, preschedule = TRUE) {
  # in RStudio, open file descriptors surpass the limit, since each child process
  # also clones RStudio's utilities; instead abandon parallelism
  if (Sys.getenv("RSTUDIO") == "1") {
    message("Detected call in RStudio, reverting to sequential execution.")
    cores <- 1
  }
  result <- parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = preschedule)
  errors <- Filter(function(x) inherits(x, "try-error"), result)
  if (length(errors) > 0) {
    stop("Parallel execution failed in ", length(errors), " worker(s): ",
         paste(vapply(errors, conditionMessage, character(1)), collapse = "; "))
  }
  result
}

# instead of modifyList() which ignores parameters with NULL values
modifyArgs <- function(original, overwritingArgs) {
  if (length(overwritingArgs) == 0) return(original)
  for (name in names(overwritingArgs)) {
    val <- overwritingArgs[[name]]
    # recurse when both are lists (and the new value is a list)
    if (is.list(original[[name]]) && is.list(val)) {
      original[[name]] <- modifyArgs(original[[name]], val)
    } else {
      # use single-bracket assignment with a list to preserve NULL as a value
      original[name] <- list(val)
    }
  }
  original
}

# Find valid ticks for a given number of bins and bp length
# Valid ticks are those where:
# 1. (ticks - 1) divides num_bins evenly
# 2. bp_length / ((num_bins / (ticks - 1)) * bin_length) is an integer
findValidTicks <- function(num_bins, bin_length, bp_length) {
  valid_ticks <- c()
  for (i in 1:num_bins) {
    if (num_bins %% i == 0) {
      # i is a valid divisor of num_bins, so ticks = i + 1
      ticks <- i + 1
      # Check if this also divides bp_length evenly
      tick_distance_bp <- (num_bins / i) * bin_length
      if (bp_length %% tick_distance_bp == 0) {
        valid_ticks <- c(valid_ticks, ticks)
      }
    }
  }
  valid_ticks
}

# Get default ticks or adjust user-provided ticks to valid ones
getValidTicks <- function(num_bins, bin_length, bp_length, user_ticks = NULL, param_name = "start_ticks") {
  valid_ticks <- findValidTicks(num_bins, bin_length, bp_length)

  if (length(valid_ticks) == 0) {
    stop("No valid tick values found for ", param_name, ". ",
         "Cannot divide bins (", num_bins, ") and bp_length (", bp_length, ") with bin_length (", bin_length, ") evenly.")
  }

  # If no user ticks provided, return a suitable default
  if (is.null(user_ticks)) {
    # Return the middle divisor for a balanced division
    return(valid_ticks[ceiling(length(valid_ticks) / 2)])
  }

  # Check if user ticks is valid
  if (user_ticks %in% valid_ticks) {
    return(user_ticks)
  }

  # Find the closest valid ticks
  closest_idx <- which.min(abs(valid_ticks - user_ticks))
  closest_ticks <- valid_ticks[closest_idx]

  warning("User-specified ", param_name, " (", user_ticks, ") is not valid for the given configuration. ",
          "Using closest valid value: ", closest_ticks, ". ",
          "Valid options are: ", paste(valid_ticks, collapse = ", "))

  return(closest_ticks)
}

# generate labels for the x axis of metagene plots
metagene_axis_labels <- function(config, region_start_label, region_end_label, start_ticks, end_ticks, 
                                 upstream_bp, downstream_bp, fixedLabelsStart = NULL, fixedLabelsEnd = NULL) {
  geneBodyBins <- config$bin_genome$bins
  if(is.null(geneBodyBins)) stop("[metagene_axis_label] Config has no bins value stored. ")
  if(geneBodyBins == 1) message("[metagene_axis_label] Only one gene body bin detected. region_start_label will be used for its tick label.")
  bodyMiddle <- switch(as.character(geneBodyBins), 
                       "1" = region_start_label, 
                       "2" = 50, 
                       "3" = c(33, 66), 
                       "4" = c(25, 50, 75), c(20, 40, 60, 80))
  # parse fixed bins from config
  fixedBinNumberUpstream <- 0 # # fixed bins (passed directly to plotGroup)
  fixedBinNumberDownstream <- 0
  fixedBinLengthUpstream <- 0
  fixedBinLengthDownstream <- 0
  upstreamLength <- 0
  downstreamLength <- 0
  # generate upstream labels
  if(length(config$bin_genome$fixedBinSizeUpstream) > 0) {
    fixedBinsUpstream <- strsplit(as.character(config$bin_genome$fixedBinSizeUpstream), ":")
    fixedBinLengthUpstream <- as.integer(fixedBinsUpstream[[1]][1])
    fixedBinNumberUpstream <- as.integer(fixedBinsUpstream[[1]][2])
    upstreamLength <- fixedBinLengthUpstream * fixedBinNumberUpstream
    
    start_ticks <- getValidTicks(fixedBinNumberUpstream, fixedBinLengthUpstream, 
                                 ifelse(is.null(upstream_bp), fixedBinNumberUpstream * fixedBinLengthUpstream, upstream_bp), 
                                 start_ticks, "start_ticks")
  }
  # generate downstream labels
  if(length(config$bin_genome$fixedBinSizeDownstream) > 0) {
    fixedBinsDownstream <- strsplit(as.character(config$bin_genome$fixedBinSizeDownstream), ":")
    fixedBinLengthDownstream <- as.integer(fixedBinsDownstream[[1]][1])
    fixedBinNumberDownstream <- as.integer(fixedBinsDownstream[[1]][2])
    downstreamLength <- fixedBinLengthDownstream * fixedBinNumberDownstream
    
    end_ticks <- getValidTicks(fixedBinNumberDownstream, fixedBinLengthDownstream, 
                                 ifelse(is.null(downstream_bp), fixedBinNumberDownstream * fixedBinLengthDownstream, downstream_bp), 
                                 end_ticks, "end_ticks")
  }
  firstBodyLabel <- "0%"
  endBodyLabel <- "100%"
  
  # validate up and downstream_bp (part of fixed regions where TSS/TTS labels are placed)
  if(fixedBinNumberUpstream > 0 && !is.null(upstream_bp) && upstream_bp > upstreamLength) 
    stop("upstream_bp must be <= fixed upstream bins * size from binning")
  if(fixedBinNumberDownstream > 0 && !is.null(downstream_bp) && downstream_bp > downstreamLength) 
    stop("downstream_bp must be <= fixed downstream bins * size from binning")
  
  # autogenerate fixed up/sownstream labels if not set
  # upstream labels
  if (fixedBinNumberUpstream > 0 && is.null(fixedLabelsStart)) {
    # if upstream_bp not given use info from fixedUpstreamBins
    effective_upstream_bp <- if (!is.null(upstream_bp)) upstream_bp else upstreamLength
    if (effective_upstream_bp < upstreamLength) {
      # calc tick dist
      startTickDistanceBins <- fixedBinNumberUpstream / (start_ticks - 1)
      startTickDistanceBp <- startTickDistanceBins * fixedBinLengthUpstream
      startTicksBeforeRegion <- as.integer(effective_upstream_bp / startTickDistanceBp)
      startTicksInRegion <- start_ticks - startTicksBeforeRegion - 1
      
      fixedLabelsStart <- c(
        seq(-effective_upstream_bp, -startTickDistanceBp, length.out = startTicksBeforeRegion),
        region_start_label,
        seq(startTickDistanceBp, upstreamLength - effective_upstream_bp, length.out = startTicksInRegion)
      )
    } else { # not extend into gene body
      firstBodyLabel <- ""
      upstream_kb <- effective_upstream_bp / 1000
      up_by <- upstream_kb / (start_ticks - 1)
      fixedLabelsStart <- c(
        paste0(seq(-upstream_kb, -up_by, by = up_by), "kb"),
        region_start_label)
    }
  }
  # downstream labels
  if (fixedBinNumberDownstream > 0 && is.null(fixedLabelsEnd)) {
    # fall back to full binned downstream length when downstream_bp not supplied
    effective_downstream_bp <- if (!is.null(downstream_bp)) downstream_bp else downstreamLength
    if (effective_downstream_bp < downstreamLength) {
      endTickDistanceBins <- fixedBinNumberDownstream / (end_ticks - 1)
      endTickDistanceBp <- endTickDistanceBins * fixedBinLengthDownstream
      endTicksAfterRegion <- as.integer(effective_downstream_bp / endTickDistanceBp)
      endTicksInRegion <- end_ticks - endTicksAfterRegion - 1
      
      fixedLabelsEnd <- c(
        seq(- (downstreamLength - effective_downstream_bp), -endTickDistanceBp, length.out = endTicksInRegion),
        region_end_label,
        seq(endTickDistanceBp, effective_downstream_bp, length.out = endTicksAfterRegion)
      )
    } else {
      endBodyLabel <- ""
      downstream_kb <- effective_downstream_bp / 1000
      down_by <- downstream_kb / (end_ticks - 1)
      fixedLabelsEnd <- c(region_end_label,
        paste0(seq(down_by, downstream_kb, by = down_by), "kb"))
    }
  }
  
  # body labels
  if(geneBodyBins == 1) {
    bodyLabels <- "" # will be ignored but has to be set
  } else {
    bodyLabels <- c(firstBodyLabel, paste0(bodyMiddle, "%"), endBodyLabel)
  }
  
  list(labels = bodyLabels, fixedLabelsStart = fixedLabelsStart, fixedLabelsStartTotalBins = fixedBinNumberUpstream, 
       fixedLabelsEnd = fixedLabelsEnd, fixedLabelsEndTotalBins = fixedBinNumberDownstream)
}

getFuncName <- function(aggr_rep_FUN) {
  if(identical(aggr_rep_FUN, base::mean)) return("mean")
  if(identical(aggr_rep_FUN, stats::median)) return("median")
  paste0("func_", substr(digest::digest(aggr_rep_FUN), 1, 7))
}



#' Heatmap of binMatrix conditions
#' 
#' @param output_pdf Path to the output PDF file.
#' @param metagene_data Output of \code{\link{load_metagene_data}}. Provide
#'   either this or \code{run_dir}.
#' @param run_dir The run folder where binning took place. Forwarded to
#'   \code{\link{load_metagene_data}} together with \code{...}.
#' @param aggr_rep_FUN Function to aggregate replicates. Default \code{mean}.
#' @param ids_to_subset Character vector of IDs to filter regions.
#' @param gsub_name Character string to remove from filenames for cleaner names.
#' @param threads Number of threads for parallel processing.
#' @param overwrite_rds_objects Logical. Overwrite cached RDS? Default \code{FALSE}.
#' @param cov_files_to_ignore Regex to exclude coverage files. Default \code{".norm.coverage.csv"}.
#' 
#' @param row_order_by How to order rows. One of:
#'   \describe{
#'     \item{\code{"mean"}}{Sort descending by row mean. Fast, O(n log n).
#'       Recommended for large matrices (>5k rows).}
#'     \item{\code{"hclust"}}{Hierarchical clustering (Ward's method) via
#'       \code{fastcluster::hclust()} + \code{parallelDist::parDist()}.
#'       Slow for large matrices (O(n²) memory and time).}
#'     \item{\code{"none"}}{Keep original row order.}
#'   }
#'   Default \code{"mean"}.
#' @param cluster_cores Integer. Number of cores for parallel distance matrix
#'   computation via \code{parallelDist::parDist()}. Only used when
#'   \code{row_order_by = "hclust"}. Default \code{1}.
#' @param id_sets Named list of character vectors containing feature ids (e.g. from \code{\link{cluster_metagene}$clusters}). 
#'  Heatmaps are annotated with these clusters, features not in one of these sets are dropped. 
#'  With \code{id_sets}, rows are ordered by \code{row_order_by} for the first condition, the order is kept for all following conditions. Set \code{"none"} for rows to be ordered exactly as the vectors.
#'  Default \code{NULL}: no subsets used, rows are ordered individually for each condition.
#' @param colors Character vector of at least 2 colors defining the gradient.
#'   Default is a viridis palette (dark purple → yellow), which maps
#'   near-zero values to dark purple rather than white/light colors.
#' @param quantile_range Numeric vector of length 2 in \code{[0, 1]} giving the
#'   lower and upper quantiles used to anchor the color scale. Values outside
#'   this range are clipped to the nearest color. Default \code{c(0.01, 0.99)},
#'   which prevents a few high-coverage outliers from washing out the rest.
#'   Set to \code{c(0, 1)} to use the full min–max range.
#' @param upstream_label,downstream_label Labels for the TSS and TTS boundary
#'   marks drawn below the heatmap. Positions are parsed automatically from the
#'   \code{@@binLength} slot. If the slot does not encode upstream/downstream
#'   regions no marks are added. Defaults: \code{"TSS"}, \code{"TTS"}.
#' @param legend_name Title of the color legend. Default \code{NULL} derives it from the
#'  normalization state recored in the binMatrix object, similar to the metagene plots.
#' @param pdf_width,pdf_height PDF dimensions in inches. Defaults: 7 x 10.
#' @param condition_mapper Forwarded to \code{\link{load_metagene_data}} when \code{run_dir} is given.
#'
#' @details
#' One page per condition. Rendering uses \code{ComplexHeatmap} with
#' \code{use_raster = TRUE} and \code{raster_resize_mat = TRUE}. Resolution is reduced 
#' to keep the PDF small for genome-wide matrices. 
#' Color scale is shared across all conditions. Row ordering is applied before
#' passing to \code{ComplexHeatmap} with \code{cluster_rows = FALSE}.
#' Column sections (upstream / gene body / downstream) are derived from the
#' \code{@@binLength} slot and annotated via \code{anno_mark}.
#' With \code{id_sets} rows are grouped into clusters in the order of the list.
#'
#' @returns Invisibly returns \code{output_pdf}.
#' @export
#'
#' @examples
#' \dontrun{
#'   heatmap_binMatrix(
#'     run_dir      = "run/",
#'     output_pdf   = "plots/heatmaps.pdf",
#'     row_order_by = "mean"
#'   )
#' }
heatmap_binMatrix <- function(output_pdf, run_dir = NULL, metagene_data = NULL,
                              condition_mapper = NULL, aggr_rep_FUN = mean, ids_to_subset = NULL, 
                              gsub_name = NULL, threads = 1, overwrite_rds_objects = FALSE, 
                              cov_files_to_ignore = ".norm.coverage.csv",
                              row_order_by = c("mean", "hclust", "none"),
                              cluster_cores = 1L,
                              id_sets = NULL,
                              colors = NULL,
                              quantile_range = c(0.01, 0.99),
                              upstream_label = "TSS",
                              downstream_label = "TTS",
                              legend_name = NULL,
                              pdf_width = 7,
                              pdf_height = 10) {
  row_order_by <- match.arg(row_order_by)
  
  # check arguments
  if(is.null(run_dir) && is.null(metagene_data)) stop("Provide either run_dir or metagene_data.")
  if(!is.null(run_dir) && !is.null(metagene_data)) stop("Provide either run_dir or metagene_data, not both.")

  # load data
  if(is.null(metagene_data)) {
    metagene_data <- load_metagene_data(
      run_dir = run_dir, condition_mapper = condition_mapper, aggr_rep_FUN = aggr_rep_FUN,
      ids_to_subset = ids_to_subset, gsub_name = gsub_name, threads = threads,
      overwrite_rds_objects = overwrite_rds_objects, cov_files_to_ignore = cov_files_to_ignore
    )
  }
  
  # validate id sets if set
  if(!is.null(id_sets)) {
    if (!is.list(id_sets) || length(id_sets) == 0 || is.null(names(id_sets)) || any(names(id_sets) == "") ||
        anyDuplicated(names(id_sets)) || !all(vapply(id_sets, is.character, logical(1))))
      stop("<id_sets> must be a list of character vectors with feature ids and unique names like $clusters from cluster_metagene()")
    if(anyNA(unlist(id_sets))) stop("<id_sets> contain NA")
    if(anyDuplicated(unlist(lapply(id_sets, unique)))) stop("<id_sets> overlap, a feature can only be in one list")
  }
  
  grouped <- metagene_data
  if (!is.list(grouped) || is.null(names(grouped))) stop("metagene_data must be a named list of binMatrix objects (see load_metagene_data()).")
  n <- length(grouped)
  if(n == 0) stop("metagene_data is empty.")
  message("[heatmap_binmatrix] Data has ", n, " condition(s).")
  
  # render title
  if(is.null(legend_name)) {
    binMat1 <- grouped[[1]]
    legend_name <- if(isTRUE(getAttribute(binMat1, "norm.binlength"))) "coverage per bp" else "coverage per bin"
    if(isTRUE(getAttribute(binMat1, "norm.libsize"))) legend_name <- paste0(legend_name, "\n/ library size * \n", 10^9)
    if (isTRUE(getAttribute(binMat1, "norm.maxShape"))) legend_name <- "coverage normed\nby shape max"
    if (isTRUE(getAttribute(binMat1, "norm.sumShape"))) legend_name <- "coverage normed\nby shape sum"
  }
  
  # package checks
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' required. Install via: BiocManager::install('ComplexHeatmap')")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' required. Install via: install.packages('circlize')")
  if (row_order_by == "hclust") {
    if (!requireNamespace("fastcluster", quietly = TRUE))
      stop("Package 'fastcluster' required for hclust. Install via: install.packages('fastcluster')")
    if (!requireNamespace("parallelDist", quietly = TRUE))
      stop("Package 'parallelDist' required for hclust. Install via: install.packages('parallelDist')")
  }

  if (is.null(colors)) colors <- grDevices::hcl.colors(9, "viridis")

  mats <- lapply(grouped, function(bm) bm@cov)
  if(!is.null(id_sets)) {
    keep <- unique(unlist(id_sets))
    mats <- lapply(mats, function(m) m[rownames(m) %in% keep, , drop = FALSE])
  }

  # Shared color function across all conditions, clipped to quantile range
  all_vals <- unlist(lapply(mats, as.vector), use.names = FALSE)
  all_vals <- all_vals[is.finite(all_vals)]
  col_limits <- quantile(all_vals, probs = quantile_range)
  message("[heatmap_binMatrix] Color scale clipped to [",
          round(col_limits[1], 3), ", ", round(col_limits[2], 3),
          "] (quantiles ", quantile_range[1], "-", quantile_range[2], ")")
  col_fun  <- circlize::colorRamp2(
    seq(col_limits[1], col_limits[2], length.out = length(colors)),
    colors
  )

  # Parse upstream/downstream bin counts from @binLength of first non-empty condition
  boundary_bins <- NULL
  for (bm in grouped) {
    bl <- bm@binLength
    if (length(bl) == 0L) next
    parts <- strsplit(bl[1L], ",")[[1L]]
    if (length(parts) < 3L) break
    up_bins   <- suppressWarnings(as.integer(strsplit(parts[1L], ":")[[1L]][2L]))
    down_bins <- suppressWarnings(as.integer(strsplit(parts[3L], ":")[[1L]][2L]))
    if (!anyNA(c(up_bins, down_bins))) {
      body_bins <- ncol(bm@cov) - up_bins - down_bins
      boundary_bins <- list(upstream = up_bins, body = body_bins, downstream = down_bins)
      message("[heatmap_binMatrix] Region boundaries: ",
              up_bins, " upstream / ", body_bins, " body / ", down_bins, " downstream bins.")
    }
    break
  }
  
  # generate cluster bands
  row_split <- NULL
  if(!is.null(id_sets)) {
    ref <- mats[[1]]
    bands <- lapply(id_sets, function(ids) rownames(order_heatmap_rows(ref[intersect(ids, rownames(ref)), , drop = FALSE], row_order_by, cluster_cores)))
    sizes <- vapply(bands, length, integer(1))
    requested <- vapply(id_sets, function(ids) length(unique(ids)), integer(1))
    if (any(sizes == 0)) stop("[heatmap_binMatrix] None of the ids of <id_sets> '", paste(names(id_sets)[sizes == 0], collapse = "', '"), "' are in the data, check format.")
    message("[heatmap_binMatrix] Heatmap bands (ids found / requested): ", paste0(names(id_sets), " = ", sizes, "/", requested, collapse = ", "),
            "; rows ordered by ", row_order_by, " of ", names(mats)[1])
    band_labels <- paste0(names(id_sets), "\n(n=", sizes, ")")
    row_split <- factor(rep(band_labels, sizes), levels = band_labels)
    row_order <- unlist(bands, use.names = FALSE)
  }

  dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)

  pdf(output_pdf, width = pdf_width, height = pdf_height)
  on.exit(dev.off(), add = TRUE)

  for (cond_name in names(mats)) {
    message("[heatmap_binMatrix] Processing condition: ", cond_name)
    mat <- mats[[cond_name]]

    if (nrow(mat) == 0L) {
      warning("Condition '", cond_name, "': empty matrix, skipping.")
      next
    }

    if(is.null(id_sets)) {
      message("[heatmap_binMatrix] ordering ", nrow(mat), " rows by ", row_order_by, " for: ", cond_name)
      mat <- order_heatmap_rows(mat, row_order_by, cluster_cores)
    } else {
      mat <- mat[row_order, , drop = FALSE]
    }

    boundary_at     <- c(boundary_bins$upstream, boundary_bins$upstream + boundary_bins$body)
    boundary_labels <- c(upstream_label, downstream_label)

    top_anno <- if (!is.null(boundary_bins))
      ComplexHeatmap::HeatmapAnnotation(
        which = "column",
        mark  = ComplexHeatmap::anno_mark(
          at        = boundary_at,
          labels    = boundary_labels,
          side      = "top",
          labels_gp = grid::gpar(fontsize = 9)
        )
      )
    else NULL

    bottom_anno <- if (!is.null(boundary_bins))
      ComplexHeatmap::HeatmapAnnotation(
        which = "column",
        mark  = ComplexHeatmap::anno_mark(
          at        = boundary_at,
          labels    = boundary_labels,
          side      = "bottom",
          labels_gp = grid::gpar(fontsize = 9)
        )
      )
    else NULL

    ht <- ComplexHeatmap::Heatmap(
      mat,
      col                = col_fun,
      na_col             = grDevices::hcl.colors(1, "viridis"),
      name               = legend_name,
      column_title       = cond_name,
      show_row_names     = FALSE,
      show_column_names  = FALSE,
      cluster_rows       = FALSE,
      cluster_columns    = FALSE,
      row_split          = row_split,
      cluster_row_slices = FALSE,
      row_title_rot      = 0,
      row_gap            = grid::unit(2, "mm"),
      top_annotation     = top_anno,
      bottom_annotation  = bottom_anno,
      use_raster         = TRUE,
      raster_resize_mat  = TRUE,
      height             = grid::unit(pdf_height - 2, "inches")
    )

    ComplexHeatmap::draw(ht)
    message("[heatmap_binMatrix] Page written for: ", cond_name)
  }

  message("[heatmap_binMatrix] Done. PDF saved: ", output_pdf)
  invisible(output_pdf)
}

#' Compare binmatrix objects
#'
#' @rdname comparator_metrics
#'
#' @param binMatrix_before,binMatrix_after An object of class binMatrix
#'
#' @returns A 3-element named integer vector corresponding to the numbers of common
#' features (intersection), of features found only in first and second matrices,
#' respectively, with names `c("commonFeats", "onlyInBefore", "onlyInAfter")`.
#' If returnFeatures = TRUE, a list is returned with the features themselves instead
#' of the numbers.
#'
#' @keywords internal
compare_features <- function(binMatrix_before, binMatrix_after, returnFeatures = F) {
  cov1 <- binMatrix_before@cov
  cov2 <- binMatrix_after@cov

  if (ncol(cov1) != ncol(cov2)) stop("Number of bins are not equal.")

  feat_names1 <- rownames(cov1)
  feat_names2 <- rownames(cov2)

  common_names <- intersect(feat_names1, feat_names2)

  firstOnly <- setdiff(feat_names1, common_names)
  secondOnly <- setdiff(feat_names2, common_names)

  if (returnFeatures) return(list(commonFeats = common_names,
                               onlyInBefore = firstOnly,
                               onlyInAfter = secondOnly))
  else return(c(commonFeats = length(common_names),
           onlyInBefore = length(firstOnly),
           onlyInAfter = length(secondOnly)))
}


#' Compare binmatrix objects
#'
#' @rdname comparator_metrics
#'
#' @param binMatrix_before An object of class binMatrix
#' @param binMatrix_after An object of class binMatrix
#' @param with_direction A logical: Should the direction of the differences
#' also be considered? If TRUE the signs are assigned from the absolute differences:
#' `y_binMatrix_before - y_binMatrix_after`.
#'
#' @returns **sq_y_dist_bins** | A numeric vector of bin-wise differences in the y values
#' of the binmatrices as in their metagene plots.
#'
#' @keywords internal
#'
sq_y_dist_bins <- function(binMatrix_before, binMatrix_after, with_direction = FALSE) {
  if (ncol(binMatrix_before@cov) != ncol(binMatrix_after@cov)) stop("Number of bins are not equal.")

  y1 <- getSummedShape(binMatrix_before)@attributes$meta
  y2 <- getSummedShape(binMatrix_after)@attributes$meta

  sign_factor <- 1
  difference <- y2-y1
  if (with_direction) sign_factor <- sign(difference)
  return(sign_factor * difference^2)
}

#' @rdname comparator_metrics
#'
#' @returns **mean_sqdist_commonFeatures** | A numeric vector of length 1.
#' Mean of the squared distances between the two binmatrices.
#' @keywords internal
mean_sqdist_commonFeatures <- function(binMatrix_before, binMatrix_after) {
  return(mean(sqdist_commonFeatures(binMatrix_before, binMatrix_after)))
}

#' @rdname comparator_metrics
#'
#' @returns **wassersteinOverBins** | A numeric vector of same length as bin number.
#' @keywords internal
wassersteinOverBins <- function(binMatrix_before, binMatrix_after) {
  cov1 <- binMatrix_before@cov
  cov2 <- binMatrix_after@cov

  if (ncol(cov1) != ncol(cov2)) stop("Number of bins are not equal.")

  dist <- sapply(1:ncol(cov1), function(col_ind) {
    wasserstein(cov1[, col_ind], cov2[, col_ind])
  })
  names(dist) <- colnames(cov1)
  return(dist)
}

#' compare metagene rpofiles of two conditions
#' 
#' This method quantifies per-bin differences between two conditions of an experiment (sample-vs-sample comparisons, for the effect of filtering on a single sample see \code{\link{filtering_impact}})
#' condition-level data returned by \code{\link{load_metagene_data}} is always normalized by library size and binlength so values are comparable across conditions by design
#' result magnitues are data-dependent (compare bins within one binning run rather than across different experiments)
#' @section Choosing a method:
#' \describe{
#'   \item{\code{"difference"}}{Signed per-bin difference of the aggregated
#'     curves (\code{conditions[2] - conditions[1]}), in the units that
#'     \code{\link{plot_metagene_experiment}} draws. Use to locate where the
#'     plotted curves diverge.}
#'   \item{\code{"rmsd"}}{Per-bin root-mean-square difference across features
#'     common to both conditions, paired by feature ID. Use when both
#'     conditions cover the same features and you want the magnitude of
#'     per-feature change.}
#'   \item{\code{"wasserstein"}}{Per-bin Wasserstein (earth mover's) distance
#'     between the distributions of coverage across features. Ignores feature
#'     identity, so it also works when the feature sets differ. Units are
#'     coverage units (the average per-feature shift needed to turn one
#'     distribution into the other); values depend on normalization and
#'     replicate count.}
#'   \item{\code{"wilcoxon"}}{Per-bin paired Wilcoxon signed-rank test across
#'     common features - the same test as the p-value strips of
#'     \code{\link{plot_metagene_experiment}}. \code{value} holds the median
#'     paired difference (effect size) and \code{n_pairs} the number of
#'     non-tied pairs actually informing the test.}
#' }
#'
#' @param metagene_data Output of \code{\link{load_metagene_data}}. Provide
#'   either this or \code{run_dir}.
#' @param run_dir The run folder where binning took place. Forwarded to
#'   \code{\link{load_metagene_data}} together with \code{...}.
#' @param conditions Character vector of length 2 naming the two conditions to
#'   compare. May be omitted when the data holds exactly two conditions.
#' @param method Comparison method, see section below. Default \code{"difference"}.
#' @param normByShapeSum,normByShapeMax Logical. Normalize each feature's shape
#'   before comparing. Defaults match \code{\link{plot_metagene_experiment}}.
#' @param aggregateFun Function aggregating features per bin
#'   (\code{method = "difference"} only). Default \code{mean}.
#' @param p_adjust Multiple-testing correction for \code{method = "wilcoxon"};
#'   any method accepted by \code{\link[stats]{p.adjust}}. Default
#'   \code{"bonferroni"} (matching the plot strips); \code{"BH"} is a less
#'   conservative alternative.
#' @param ... Forwarded to \code{load_metagene_data} when \code{run_dir} is used.
#'
#' @returns A data.frame with one row per bin and a fixed schema: \code{bin},
#'   \code{cond_a}, \code{cond_b}, \code{method}, \code{value}, \code{p_raw},
#'   \code{p_adj}, \code{n_pairs}. The last three are \code{NA} except for
#'   \code{method = "wilcoxon"}.
#' 
#' @export
compare_metagene <- function(metagene_data = NULL, run_dir = NULL, conditions = NULL, 
                             method = c("difference", "rmsd", "wasserstein", "wilcoxon"), 
                             normByShapeSum = FALSE, normByShapeMax = FALSE, 
                             aggregateFun = mean, p_adjust = "bonferroni", ...) {
  method <- match.arg(method)
  if(is.null(run_dir) && is.null(metagene_data)) stop("Provide either run_dir or metagene_data")
  if(!is.null(run_dir) && !is.null(metagene_data)) stop("Provide either run_dir or metagene_data, not both.")
  if(!is.null(run_dir)) {
    metagene_data <- load_metagene_data(run_dir = run_dir, ...)
  }
  if(is.null(conditions)) {
    if(length(metagene_data) != 2) stop("can only compare two conditions, specify them via <conditions>. Available conditions: ", paste(names(metagene_data), collapse = ", "))
    conditions <- names(metagene_data)
  }
  if(length(conditions) != 2) stop("conditions must name exactly two conditions")
  missing_conds <- setdiff(conditions, names(metagene_data))
  if(length(missing_conds)>0) stop("Conditions not in metagene data: ", paste(missing_conds, collapse = ", "), 
                                   ". Available conditions: ", paste(names(metagene_data), collapse = ", "))
  message("[compare_metagene] ", conditions[1], " vs. ", conditions[2], ", method: ", method)
  shape_1 <- getSummedShape(metagene_data[[conditions[1]]], normByShapeSum = normByShapeSum, normByShapeMax = normByShapeMax, aggregateFun = aggregateFun)
  shape_2 <- getSummedShape(metagene_data[[conditions[2]]], normByShapeSum = normByShapeSum, normByShapeMax = normByShapeMax, aggregateFun = aggregateFun)
  cov_1 <- shape_1@attributes[["cov"]]
  cov_2 <- shape_2@attributes[["cov"]]
  if(ncol(cov_1) != ncol(cov_2)) stop("bin numbers differ between conditions (", ncol(cov_1), " vs ", ncol(cov_2), ").")
  nr_bins <- ncol(cov_1)
  res <- data.frame(bin = seq_len(nr_bins), cond1 = conditions[1], cond2 = conditions[2], method = method)
  if(method == "difference") {
    res$value <- shape_2@attributes[["meta"]] - shape_1@attributes[["meta"]]
  } else if(method == "rmsd") {
    paired <- pair_by_rownames(cov_1, cov_2)
    res$value <- sqrt(colMeans((paired$b - paired$a)^2))
  } else if(method == "wasserstein") {
    res$value <- vapply(seq_len(nr_bins), function(i) wasserstein(cov_1[, i], cov_2[, i]), numeric(1))
  } else if(method == "wilcoxon") {
    wilcox_res <- binwise_wilcox(cov_1, cov_2, p_adjust = p_adjust)
    res$value <- wilcox_res$effect
    res$p <- wilcox_res$praw
    res$padj <- wilcox_res$padj
    res$npairs <- wilcox_res$npairs
  }
  return(res)
}


#' Quantify the impact of filtering on a metagene profile
#'
#' Diagnostic report for the before/after comparison of one sample: how
#' strongly (and where) did a processing step change the aggregated metagene
#' curve, and what did the removed features look like? For comparing two
#' different samples see \code{\link{compare_metagene}}.
#'
#' The function detects what kind of change it is looking at:
#' \describe{
#'   \item{\code{"filtering"}}{\code{after} features are a subset of
#'     \code{before} and the kept rows are unchanged. Full report, including
#'     the removed-set profile and the share of signal removed.}
#'   \item{\code{"values_changed"}}{Identical features, different values
#'     (an upstream processing change). Percent-change report only.}
#'   \item{\code{"population_changed"}}{Feature sets differ beyond a subset
#'     relation (e.g. gene-level vs transcript-level annotation). Curve-level
#'     comparison only, since changes cannot be attributed to individual rows.}
#' }
#'
#' To compare the effect of a plotting normalization instead, call
#' \code{\link{get_metagene_profiles}} twice with different norm settings.
#'
#' @param before,after Outputs of \code{\link{load_metagene_data}} holding the
#'   same condition before and after the processing step.
#' @param condition Name of the condition to compare. May be omitted when both
#'   objects hold exactly one condition.
#' @param aggregateFun Function aggregating features per bin. Default \code{mean}.
#' @returns A list with \code{summary} (comparison_type, feature counts,
#'   percent of features and of signal removed, maximum absolute percent
#'   change and its bin, number of zero-baseline bins) and \code{per_bin}
#'   (data.frame: \code{bin}, \code{curve_before}, \code{curve_after},
#'   \code{pct_change}, \code{curve_removed}).
#' @export
filtering_impact <- function(before, after, condition = NULL, aggregateFun = mean) {
  # check validity of inputs
  # only one condition present or <condition> set? is the condition present in both before/after?
  if(is.null(condition)) {
    if(length(before) != 1) stop("<before> has more than one condition, either subset <before> or provide <condition>")
    if(length(after) != 1) stop("<after> has more than one condition, either subset <after> or provide <condition>")
    binmat_before <- before[[1]]
    binmat_after <- after[[1]]
  } else {
    if(!condition %in% names(before)) stop("<condition> not found in <before>. Available conditions: ", paste(names(before), collapse = ", "))
    if(!condition %in% names(after)) stop("<condition> not found in <after>. Available conditions: ", paste(names(after), collapse = ", "))
    binmat_before <- before[[condition]]
    binmat_after <- after[[condition]]
  }
  # check bin layout, normalization, duplicates
  if(ncol(binmat_before@cov) != ncol(binmat_after@cov)) stop("Matrices are not comparable, bin numbers differ!")
  norm_attributes <- c("norm.libsize", "norm.binlength", "norm.sumShape", "norm.maxShape")
  if(!identical(binmat_before@attributes[norm_attributes], binmat_after@attributes[norm_attributes])) 
    stop("Matrices are not comparable, normalization differs!")
  rows_before <- rownames(binmat_before@cov)
  rows_after <- rownames(binmat_after@cov)
  if(anyDuplicated(rows_before) || anyDuplicated(rows_after)) stop("duplicate rownames detected, comparison would be ambiguous")
  
  # which checks should be done? filtering -> full report, values_changed -> report percentage change, population_changed -> (eg gene vs transcript) curve level comparison
  common_rows <- intersect(rows_before, rows_after)
  is_unchanged <- length(common_rows) > 0 && isTRUE(all.equal(binmat_before@cov[common_rows,], binmat_after@cov[common_rows,]))
  if(setequal(rows_before, rows_after)) {
    if(is_unchanged) {
      comparison_type <- "identical"
      warning("<before> and <after> are identical, nothing changed")
    } else {
      comparison_type <- "values_changed"
    }
  } else if(all(rows_after %in% rows_before) && is_unchanged) {
    comparison_type <- "filtering"
  } else {
    comparison_type <- "population_changed"
    message("Feature sets differ (", round(100*length(common_rows)/max(1, length(rows_before)), 1), "% overlap)")
  }
  
  curve_before <- apply(binmat_before@cov, 2, aggregateFun)
  curve_after <- apply(binmat_after@cov, 2, aggregateFun)
  percent_changed <- ifelse(curve_before == 0, NA, 100*(curve_after - curve_before) / curve_before)
  
  removed <- setdiff(rows_before, rows_after)
  curve_removed <- rep(NA, length(curve_before))
  percent_signal_removed <- NA
  if(comparison_type == "filtering") {
    curve_removed <- apply(binmat_before@cov[removed, , drop = FALSE], 2, aggregateFun)
    percent_signal_removed <- 100* sum(binmat_before@cov[removed, , drop = FALSE])/(sum(binmat_before@cov))
  }
  max_change_bin <- ifelse(all(is.na(percent_changed)), NA, which.max(abs(percent_changed)))
  
  results <- list(summary = list(
    comparison_type = comparison_type, 
    n_before = length(rows_before), 
    n_after = length(rows_after), 
    n_removed = length(removed), 
    percent_features_removed = 100 * length(removed) / length(rows_before), 
    percent_signal_removed = percent_signal_removed, 
    largest_percent_change = ifelse(is.na(max_change_bin), NA, percent_changed[max_change_bin]), 
    largest_percent_change_atBin = max_change_bin
  ), 
  per_bin = data.frame(bin = seq_along(curve_before), 
                       curve_before = curve_before, 
                       curve_after = curve_after, 
                       percent_changed = percent_changed, 
                       curve_removed = curve_removed))
  return(results)
}

#' Detect outlier ids across metagene conditions
#' 
#' Flags features that have extreme coverage by three different criteria: 
#' \describe{
#'  \item{\code{influence}}{ (default) Use this to check whether peaks in the metagene profile are driven by single features. 
#'  The total coverage is split among all features at each bin. Influence of a feature is then defined as 
#'  the largest percentage of coverage the feature has in any bin. \code{threshold = 0.1} flags features that account for more than 10% 
#'  of the signal at some point of the curve, while the average feature accounts for \code{1/n} as printed in the report for reference. }
#'  \item{\code{MAD}}{Use this to check whether a feature has extreme coverage compared to average features in this set. 
#'  Reports a features maximum bin coverage and flags if it is more than \code{threshold = 3} MADs above the median. 
#'  Note: coverage maxima can be right-skewed, increase threshold after inspecting the report. }
#'  \item{\code{quantile}}{Use this to rank cut top features. 
#'  Features are ordered by their maximum bin coverage and flagged if they are above the \code{threshold = 0.99} quantile. }
#' }
#' Returns ID set \code{$keep} for passing \code{ids_to_subset} to \code{load_metagene_data()}, 
#' \code{plot_metagene_experiment()}, \code{get_metagene_profiles()} or \code{plot_metagene_from_json()}
#' to apply the filtering. 
#' 
#' @param binMat A single \code{binMatrix} or a named list of \code{binMatrix} objects as generated from \code{load_metagene_data()}. 
#'  For a list, outlier detection runs per condition and the union of IDs is flagged for removal
#' @param method \code{influence}, \code{MAD}, \code{quantile}
#' @param threshold Numeric scalar with specific meaning per method (see above). 
#'  Defaults: 0.1 for influence (allowed values ]0, 1[), 3 for MAD (allowed values > 0), 0.99 for quantile (allowed values ]0, 1[). 
#'  Pass \code{NA} for report-only mode (nothing is flagged, \code{stats} can be used to choose a threshold).
#' 
#' @returns A list with \code{keep} and \code{remove} (two character vecotrs with feature IDs, remove is the union of flagged features over conditions),
#'  \code{stats} (named list with one named numeric vector of per-feature statistics per condition), 
#'  \code{cutoffs} (named numeric vector with one value per condition), 
#'  \code{method} and \code{threshold}
#' @export
detect_outlier_features <- function(binMat, method = c("influence", "MAD", "quantile"), threshold = NULL){
  method <- match.arg(method)
  
  report_only <- FALSE
  if(is.null(threshold)) {
    threshold <- switch(method, influence = 0.1, MAD = 3, quantile = 0.99)
  } else if(is.na(threshold)){
    report_only <- TRUE
  } else if(!is.numeric(threshold) || length(threshold) != 1) {
    stop("<threshold> must be a single number or NA for report-only mode.")
  } else if(method == "MAD"){
    if(threshold <= 0) stop("<threshold> must be positive.")
  } else if(threshold <= 0 || threshold >= 1) {
    stop("<threshold> must be in ]0, 1[")
  }
  
  if(inherits(binMat, "binMatrix")) binMat <- list(binMatrix = binMat)
  if(!is.list(binMat) || length(binMat) == 0 || !all(vapply(binMat, inherits, logical(1), what = "binMatrix")))
    stop("<binMat> must be a binMatrix or a list of binMatrix objects e.g. generated by load_metagene_data()")
  if(is.null(names(binMat))) names(binMat) <- paste0("condition_", seq_along(binMat))
  
  cutoffs <- numeric(0)
  stats_list <- list()
  removed <- character(0)
  all_ids <- character(0)
  
  #loop over binmatrices, find values
  for(cond in names(binMat)) {
    # checks
    mat <- binMat[[cond]]@cov
    if(nrow(mat) == 0) stop("BinMatrix for ", cond, " has no rows.")
    if(is.null(rownames(mat))) stop("Coverage matrix for ", cond, " has no rownames.")
    all_ids <- union(all_ids, rownames(mat))
    peak_bin <- NULL
    
    if(method == "influence"){
      if(any(!is.finite(mat))) stop("Coverage matrix for condition ", cond, " contains non-finite values.")
      col_total <- colSums(mat)
      use <- col_total > 0
      if(!any(use)) stop("All bins have 0 total coverage in condition ", cond, ".")
      fractions <- sweep(mat[, use, drop = FALSE], 2, col_total[use], "/") # what is the share of each feature in each bin
      values <- apply(fractions, 1, max)
      peak_bin <- colnames(fractions)[apply(fractions, 1, which.max)]
      cutoff <- if(report_only) NA else threshold
    } else {
      values <- apply(mat, 1, function(val) {
        val <- val[is.finite(val)]
        if(length(val) == 0) NA else max(val)
      })
      if(anyNA(values)) stop(sum(is.na(values)), " rows are NA, remove these first.")
      if(report_only) {
        cutoff <- NA
      } else if(method == "MAD") {
        median_val <- stats::median(values)
        mad_val <- stats::mad(values)
        if(mad_val <= .Machine$double.eps * max(1, abs(median_val))) # check if it is zero but with a safety net!
          stop("MAD of the per-feature maxima is zero, more than half of the features share the same value. ", 
               "Use other method or remove unexpressed features first.")
        cutoff <- median_val + threshold * mad_val
      } else {
        cutoff <- stats::quantile(values, probs = threshold, names = FALSE)
      }
    }
    
    names(values) <- rownames(mat)
    stats_list[[cond]] <- values
    cutoffs[cond] <- cutoff
    
    flagged <- if (report_only) rep(FALSE, nrow(mat)) else values > cutoff
    removed <- union(removed, names(values)[flagged])
    
    if(! report_only && sum(flagged) > 0.5 * length(flagged)) 
      warning("[detect_outlier_features] more than half of the rows for condition ", cond, " are flagged for removal.")
    
    header <- if (report_only)
      sprintf("[detect_outlier_features] %s (%s, report only): top features",
              cond, method)
    else
      sprintf("[detect_outlier_features] %s (%s, threshold = %s): %d / %d features above cutoff = %.4g",
              cond, method, threshold, sum(flagged), nrow(mat), cutoff)
    
    if(method == "influence") {
      header <- paste0(header, sprintf(" (average share 1/n = %.2g)", 1 / nrow(mat)))
    
      top10 <- order(values, decreasing = TRUE)[seq_len(min(10, nrow(mat)))] # top ten or less if matrix only has a few rows
      lines <- if(method == "influence") sprintf("  - %s (share: %.1f%% at bin %s)", names(values)[top10], 100 * values[top10], peak_bin[top10]) 
        else sprintf("  - %s (max coverage: %.4g)", names(values)[top10], values[top10])
      message(paste(c(header, lines), collapse = "\n"))
    }
  }
  kept <- setdiff(all_ids, removed)
  if(length(kept) == 0) stop("All features are flagged as outliers!")
  if(length(binMat) > 1 && !report_only) 
    message("[detect_outlier_features] union across ", length(binMat), " conditions: ", 
            length(removed), " unique ids flagged, ", length(kept), " are kept.")
  list(keep = kept, remove = removed, stats = stats_list, cutoffs = cutoffs, method = method, threshold = threshold)
}

#' cluster features by metagene profiles
#' 
#' Groups features into \code{k} clusters by the shape of their coverage profile. 
#' Each profile is scaled to a sum of 1 with \code{normBy = "sum"} per default (like normByShapeSum = TRUE for plotting) or alternatively to the shape max or not at all (only normalization to sequencing depth and binlength at load_metagene()).
#' For several conditions given in a list of binMatrices, profiles are scaled together across conditions. 
#' That is, features with a similar coverage distribution across conditions end up in the same cluster regardless of the extent of their coverage. 
#' With \code{normBy  = "none"}, coverage values are used as is and strongly covered features separate from weakly covered. 
#' 
#' Rows are clustered hierarchically with Ward's method on euclidean distances. 
#' 
#' All conditions must contain the same features in the same order, which should be given after one binning run. 
#' Features with zero coverage are dropped for clustering and returned with \code{$no_coverage}
#' 
#' @param binmat A single \code{binMatrix} or a named list of \code{binMatrix} objects as generated from \code{load_metagene_data()}
#' @param k Number of clusters, a single integer >= 2
#' @param conditions Character vector naming conditions of \code{binmat} that should be used. Default NULL
#' @param normBy Normalization applied to the profile of each feature (across all conditions)
#' \code{"sum"} default, profile sums to 1
#' \code{"max"} divides the profile by its largest value
#' \code{"none"} uses the coverage values (which are already normalized by binlength and lib size via load_metagene())
#' 
#' @returns A list with \code{$clusters}: named list with ID sets, 
#' \code{$tree}: the hclust object obtained from clustering, 
#' \code{$k}: number of clusters, 
#' \code{$no_coverage}: features with zero coverage that are dropped for clustering.
#' @export
#' 
cluster_metagene <- function(binmat, k, conditions = NULL, normBy = c("sum", "max", "none"), cluster_cores = 1){
  # get params, check them
  normBy <- match.arg(normBy)
  if(!is.numeric(k) || length(k) != 1 || !is.finite(k) || k < 2 || k != round(k)) stop("<k> must be a finite number >= 2")
  if(inherits(binmat, "binMatrix")) binmat <- list(condition_1 = binmat)
  if(!is.list(binmat) || length(binmat) == 0 || !all(vapply(binmat, inherits, logical(1), what = "binMatrix"))) 
    stop("<binmat> must be a binMatrix object or a list of binMatrix objects e.g. from load_meatgene_data()")
  if(is.null(names(binmat))) names(binmat) <- paste0("condition_", seq_along(binmat))
  
  if(!is.null(conditions)) {
    conditions <- unique(conditions)
    missing_conds <- setdiff(conditions, names(binmat))
    if(length(missing_conds) > 0) 
      stop("Conditions not in metagene data: ", paste(missing_conds, collapse = ", "), 
           ". Available conditions: ", paste(names(binmat), collapse = ", "))
    binmat <- binmat[conditions]
  }
  
  # check if packages are available
  if(!requireNamespace("fastcluster", quietly = TRUE))
    stop("Package 'fastcluster' required for clustering. Install via: install.packages('fastcluster')")
  if(!requireNamespace("parallelDist", quietly = TRUE))
    stop("Package 'parallelDist' required for clustering. Install via: install.packages('parallelDist')")
  
  # check that features and bins match across conditions
  ref <- binmat[[1]]@cov
  if(is.null(rownames(ref)) || anyDuplicated(rownames(ref))) 
    stop("coverage matrices need unique rownames (featureIDs).")
  for(cond in names(binmat)) {
    cur <- binmat[[cond]]@cov
    if(!identical(rownames(ref), rownames(cur)) || ncol(ref) != ncol(cur)) 
      stop("Binmatrices must have the same feature set and number of bins across conditions.")
    if(any(!is.finite(cur)) || any(cur < 0) ) stop("Binmatrices must contain positive, finite values.")
  }
  nr_bins <- ncol(ref)
  # gather conditions
  full <- do.call(cbind, lapply(binmat, function(binm) binm@cov))
  no_cov_rows <- rownames(full)[rowSums(full) == 0]
  if(length(no_cov_rows)) message("[cluster_metagene] Dropping ", length(no_cov_rows), " with zero coverage in all conditions.")
  full <- full[rowSums(full)>0, , drop = FALSE]
  if(length(full) <= k) stop("[cluster_metagene] Less than", k, "features left to cluster.")
  
  # normalization 
  if(normBy == "sum") full <- full / rowSums(full)
  if(normBy == "max") full <- full / apply(full, 1, max)
  
  tree <- fastcluster::hclust(parallelDist::parDist(full, method = "euclidean", threads = cluster_cores), method = "ward.D2") # faster than the heatmap thing
  tree$labels <- rownames(full)
  clusters <- stats::cutree(tree, k = k)
  cluster_list <- split(rownames(full), clusters)
  names(cluster_list) <- paste0("cluster", names(cluster_list))
  
  return(list(clusters = cluster_list, tree = tree, k = k, no_coverage = no_cov_rows))
}



################## HELPERS #####################################################

# vectors represent counts of different features in one bin;
# function quantifies the dissimilarity between count distributions,
# measuring the cost required to transform the signal intensity profile
# of one gene set into the other regardless of row number
wasserstein <- function(v1, v2) {
  if(length(v1) == 0 || length(v2) == 0) stop("cannot compute a wasserstein distance in an empty vector")
  if(anyNA(v1) || anyNA(v2)) stop("NA values were passed to wasserstein() - check the input matrices")
  
  v1 <- sort(v1)
  v2 <- sort(v2)

  n1 <- length(v1)
  n2 <- length(v2)

  # Unique values or breakpoints of both vectors
  all_vals <- sort(unique(c(v1, v2)))

  # calculate eCDF at the breakpoints
  # FindInterval returns how many elements are <= each breakpoint
  cdf1 <- findInterval(all_vals, v1) / n1
  cdf2 <- findInterval(all_vals, v2) / n2

  # calculate area between CDF curves
  diffs <- diff(all_vals)

  # Absolute difference between CDFs at each step
  heights <- abs(cdf1 - cdf2)[-length(all_vals)]

  # Sum of (width * height)
  distance <- sum(diffs * heights)

  return(distance)
}

#align two coverage matrices on common rownames
# errors on duplicates or empty intersections, warns on tiny overlap
pair_by_rownames <- function(cov_a, cov_b) {
  if(anyDuplicated(rownames(cov_a)) || anyDuplicated(rownames(cov_b))) stop("duplicate rownames detected")
  common <- intersect(rownames(cov_a), rownames(cov_b))
  if(length(common) == 0) stop("no common rownames between the two matrices, pairing needs shared features")
  if(length(common) < 10) warning("only few common features, results will be unreliable")
  list(a = cov_a[common, , drop = FALSE], b = cov_b[common, , drop = FALSE])
}

# binwise paired wilcoxon between two coverage matrices paired by rownames
# returns: praw, padj, npairs (non tied pairs) and median paired difference per bin
binwise_wilcox <- function(cov_a, cov_b, p_adjust = "bonferroni") {
  paired <- pair_by_rownames(cov_a, cov_b)
  a <- paired$a
  b <- paired$b
  nbins <- ncol(a)
  praw <- numeric(nbins)
  npairs <- integer(nbins)
  effect <- numeric(nbins)
  for(p in seq_len(nbins)) {
    d <- a[, p] - b[, p]
    npairs[p] <- sum(d != 0)
    effect[p] <- stats::median(d)
    praw[p] <- if(npairs[p] == 0) NaN else stats::wilcox.test(a[, p], b[, p], exact = FALSE, paired = TRUE)$p.val
  }
  padj <- stats::p.adjust(praw, method = p_adjust, n = nbins)
  praw[is.nan(praw)] <- 1 # like in legacy code
  padj[is.nan(padj) | is.na(padj)] <- 1
  data.frame(praw = praw, padj = padj, npairs = npairs, effect = effect)
}

order_heatmap_rows <- function(mat, row_order_by, cluster_cores){
  if(nrow(mat) < 2L || row_order_by == "none") return(mat)
  if(row_order_by == "mean") return(mat[order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE), , drop = FALSE])
  tmp <- mat
  tmp[!is.finite(tmp)] <- 0
  mat[fastcluster::hclust(parallelDist::parDist(tmp, threads = cluster_cores), method = "ward.D2")$order, , drop = FALSE]
}

#' Custom theme for plotting
#'
#' This function creates a custom theme for plotting in R. It is based on the `theme_bw` function from the `ggplot2` package,
#' with modifications to remove the panel grid lines. Additional arguments can be passed to the `theme` function to further
#' customize the plot.
#'
#' @param base_size The base font size for the theme (default is 16).
#' @param base_family The base font family for the theme (default is an empty string).
#' @param base_line_size The base line size for the theme, calculated as `base_size / 22` (default is `base_size / 22`).
#' @param base_rect_size The base rectangle size for the theme, calculated as `base_size / 22` (default is `base_size / 22`).
#' @param ... Additional arguments to be passed to the `theme` function.
#'
#' @return A modified theme object for plotting.
#'
#' @importFrom ggplot2 theme_bw theme element_blank
#'
#' @examples
#' theme_custom()
#'
#' @export
theme_custom <- function(base_size = 16,
                         base_family = "",
                         base_line_size = base_size / 22,
                         base_rect_size = base_size / 22,
                         ...) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family, base_line_size = base_size / 22, base_rect_size = base_size / 22) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...)
}


#' Embedding Plot
#'
#' Plot a dimension reduction embedding with points coloured based on a metadata variable or feature.
#'
#' @param data The data for the plot as Seurat object or data frame.
#' @param colour The variable used for colouring the plot.
#' @param ... Additional fixed aesthetics passed to the `geom_` function.
#' @param x The variable used for the x-axis of the plot. Default is "umap_1".
#' @param y The variable used for the y-axis of the plot. Default is "umap_2".
#' @param type The type of colour scale to use, options are "auto", "discrete", "sequential", and "diverging". Default is "auto".
#' @param cmap The colour map to use for the plot as a vector of colours. Optional names can be provided for the colours, which
#' will be used to assign colours to corresponding values in the `colour` variable. If not provided, a default colour map is used based
#' on the inferred data type.
#' @param limits The limits for the colour scale. Optional names can be provided for the limits, which will be used as labels if
#' `binarise` is set to `TRUE`. Percentile values for a limit can be provided as a character in the form 'pX' where 'X' is the percentile
#' to use. `NA` values for a limit will infer the limit from the data.
#' @param breaks The breaks for the colour scale. Default is NULL, which will use the default breaks for the scale.
#' @param midpoint Optional midpoint for the colour scale if `type` is "diverging". Default is 0.
#' @param binarise Whether to binarise the colour scale labels, taken from names of `limits`. Default is FALSE.
#' @param group The variable used for grouping points. If provided, points will be aggregated using the `group` variable and plotted using
#' the mean embedding coordinates. If `type` is "discrete",`colour` will be set to the most frequent value per group; if `type` is
#' "sequential" or "diverging", `colour` will be set to the mean value per group. Default is NULL.
#' @param split The variable used for faceting the plot. Default is NULL.
#' @param ncol The number of columns for the facets. If not provided, is automatically determined based on the number of facets.
#' @param size.adj The adjustment factor for point size. Default is 1.
#' @param labs A named list of arguments to pass to the `labs` function. Default is an empty list.
#' @param theme A named list of arguments to pass to the `theme` function. Default is an empty list.
#' @param rasterise Rasterise points using `geom_scattermore`, useful for rendering a large number of points. If not set, will be
#' automatically enabled when the number of points exceeds 100,000.
#' @param downsample Downsample points; if `type` is "discrete", downsampling will be performed for each unique value of `colour`.
#' Default is Inf.
#' @param shuffle Shuffle the order of points. Default is TRUE.
#' @param order Order the points by the `colour` variable. Default is the opposite of `shuffle`.
#' @param seed Random seed for reproducibility. Default is 42.
#' @param verbose Display status messages. Default is TRUE.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # Example usage with a Seurat object
#' embedding.plot(data = seurat_object, colour = "gene_expression")
#'
#' # Example usage with a data frame
#' embedding.plot(data = my_data, colour = "value", type = "sequential")
#'
#' @import ggplot2
#' @importFrom dplyr group_by slice_sample summarise ungroup n
#' @importFrom pals kovesi.linear_blue_95_50_c20 coolwarm
#' @importFrom scales squish rescale_mid
#' @importFrom scattermore geom_scattermore
#' @importFrom Seurat FetchData
#'
#' @export
embedding.plot <- function(data, colour, ..., x = "umap_1", y = "umap_2",
                           type = c("auto", "discrete", "sequential", "diverging"),
                           cmap = NULL, limits = c(Low = NA, High = NA), breaks = NULL,
                           midpoint = 0, binarise = FALSE, group = NULL, split = NULL, ncol = NULL,
                           size.adj = 1, labs = list(), theme = list(), rasterise = NULL,
                           downsample = Inf, shuffle = TRUE, order = !shuffle,
                           seed = 42, verbose = TRUE) {
  set.seed(seed)

  if (inherits(data, "Seurat")) {
    data <- Seurat::FetchData(object = data, vars = c(x, y, colour, group, split))
  }

  type <- match.arg(type)
  if (type == "auto") {
    type <- if (is.numeric(data[[colour]])) "sequential" else "discrete"
    if (verbose) message(sprintf("Inferred data type as '%s'. To override, set type explicitly.", type))
  }

  if (type == "discrete") {
    data <- data |>
      dplyr::group_by(get(colour)) |>
      dplyr::slice_sample(n = downsample) |>
      dplyr::ungroup()
  } else {
    data <- data |>
      dplyr::slice_sample(n = downsample)
    if (any(is.character(limits), na.rm = TRUE)) {
      names <- names(limits)
      limits <- sapply(
        limits,
        function(x) {
          if (!is.na(x) && is.character(x)) {
            if (!grepl(pattern = "^p", x = x)) {
              y <- NA
            } else {
              p <- as.numeric(sub(pattern = "^p", replacement = "", x = x))
              if (is.na(p) || p < 0 || p > 100) {
                y <- NA
              } else {
                y <- quantile(data[[colour]], probs = p / 100, na.rm = TRUE)
              }
            }
          } else {
            y <- x
          }
          if (!is.na(x) && is.na(y)) warning("Character limits must be in the form 'pX' where 'X' is the percentile to use. Ignoring invalid limit: ", x)
          return(as.numeric(y))
        }
      ) |> setNames(names)
    }
  }
  if (!is.null(group)) {
    suppressMessages({
      data <- if (!is.null(split)) data |> dplyr::group_by("{group}" := get(group), "{split}" := get(split)) else data |> dplyr::group_by("{group}" := get(group))
      data <- data |>
        dplyr::summarise(
          n = dplyr::n(),
          "{x}" := mean(get(x), na.rm = TRUE),
          "{y}" := mean(get(y), na.rm = TRUE),
          "{colour}" := if (type == "discrete") names(which.max(table(get(colour))))[1] else mean(get(colour), na.rm = TRUE)
        ) |>
        dplyr::ungroup()
    })
  }
  if (shuffle) data <- data[sample(nrow(data)), ]
  if (order) data <- data[order(data[[colour]]), ]
  if (shuffle && order) warning("Both shuffle and order are set to TRUE, order will take precedence.")

  if (is.null(cmap)) {
    cmap <- switch(type,
      sequential = pals::kovesi.linear_blue_95_50_c20(100),
      diverging = pals::coolwarm(100),
      NULL
    )
    if (verbose) message(sprintf("Using default colour map for type '%s'. To override, set cmap explicitly.", type))
  }

  args <- list(limits = limits, oob = scales::squish, breaks = if (binarise) function(x) x else if (is.null(breaks)) ggplot2::waiver() else breaks, labels = if (binarise) names(limits) else ggplot2::waiver())
  scale <- switch(type,
    discrete = if (is.null(cmap)) {
      ggplot2::scale_colour_hue()
    } else {
      ggplot2::scale_colour_manual(values = if (is.null(names(cmap))) cmap[seq_along(unique(data[[colour]]))] |> setNames(unique(data[[colour]])) else cmap)
    },
    sequential = do.call(ggplot2::scale_colour_gradientn, modifyList(args, list(colours = cmap), keep.null = TRUE)),
    diverging = do.call(ggplot2::scale_colour_gradientn, modifyList(args, list(colours = cmap, rescaler = function(x, to = c(0, 1), from, mid = midpoint) scales::rescale_mid(x, to, from, mid)), keep.null = TRUE))
  )

  if (verbose && missing(size.adj)) message("Auto-scaling point size. To override, use size.adj.")
  if (is.null(group)) {
    if (is.null(rasterise)) {
      rasterise <- nrow(data) > 1e5
      if (verbose && rasterise) message("Rasterising due to number of points exceeding 100,000. To override, set rasterise to FALSE.")
    }
    if (rasterise) {
      size <- size.adj # default point size of 1px for rasterised points
      args <- list(pointsize = size)
    } else {
      size <- size.adj / sqrt(nrow(data)) # heuristic for optimal point size based on number of points
      args <- list(size = size)
    }
  } else {
    rasterise <- FALSE
    size <- size.adj * 2 * round(log(max(data$n, na.rm = TRUE) - min(data$n, na.rm = TRUE))) # heuristic for optimal maximum point size based on range of n across groups
    args <- list(mapping = ggplot2::aes(size = n, shape = 16))
  }
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = get(x), y = get(y), colour = get(colour))) +
    scale +
    ggplot2::scale_size_area(max_size = size) +
    do.call(ggplot2::labs, modifyList(list(x = gsub(pattern = "_", replacement = "", x = x) |> toupper(), y = gsub(pattern = "_", replacement = "", x = y) |> toupper(), size = "Cells"), labs, keep.null = TRUE)) +
    do.call(theme_custom, modifyList(list(aspect.ratio = 1, axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()), theme, keep.null = TRUE)) +
    guides(size = ggplot2::guide_none()) +
    if (rasterise) do.call(scattermore::geom_scattermore, modifyList(args, list(...), keep.null = TRUE)) else do.call(ggplot2::geom_point, modifyList(args, list(...), keep.null = TRUE))

  if (!is.null(split)) {
    plot <- plot + ggplot2::facet_wrap(~ get(split), ncol = ifelse(is.null(ncol), ceiling(sqrt(length(unique(data[[split]])))), ncol))
  }

  return(plot)
}

#' Expression plot
#'
#' Wrapper around `embedding.plot` to show feature expression on a dimension reduction embedding with typical default settings.
#'
#' @param data The data for the plot as Seurat object or data frame.
#' @param features The name of the features to plot.
#' @param ... Additional arguments to be passed to the `embedding.plot` function.
#' @param thresholds Numeric scalar or vector of the same length as `features` specifying the lower threshold value(s)
#' for the expression data. Default is NA.
#' @param assay The name of the assay to search for features in if `data` is a Seurat object. Default is NULL
#' (if `key` is also not specified, will search default assay then search other assays if features not found).
#' @param key The key of the to search for features in if `data` is a Seurat object. Default is the key associated
#' with the specified assay.
#' @param verbose Display status messages. Default is TRUE.
#'
#' @return A list of `ggplot` objects.
#'
#' @examples
#' # Load Seurat object
#' data <- Read10X(data.dir = "path/to/data")
#' seurat_obj <- CreateSeuratObject(counts = data)
#'
#' # Plot expression of a feature
#' expression.plot(data = seurat_obj, features = "CD3E")
#'
#' @importFrom Seurat Key
#'
#' @export
expression.plot <- function(data, features, ..., thresholds = NA, assay = NULL, key = Seurat::Key(data[[assay]]), verbose = TRUE) {
  if (!inherits(data, "Seurat")) {
    if (verbose && !is.null(assay)) message("Ignoring assay and key because data is not a Seurat object.")
    assay <- key <- NULL
  } else {
    if (verbose && !is.null(assay) && !missing("key")) message("Both assay and key are provided, key will take precedence.")
  }
  plots <- setNames(nm = features) |>
    mapply(
      FUN = function(feature, threshold, key, ...) {
        do.call(
          embedding.plot,
          modifyList(
            list(
              data = data,
              colour = paste0(key, feature),
              type = "sequential",
              limits = c(Low = threshold, High = NA),
              binarise = TRUE,
              labs = list(title = feature, colour = NULL),
              verbose = verbose
            ),
            list(...),
            keep.null = TRUE
          )
        )
      },
      feature = features,
      threshold = rep_len(thresholds, length(features)),
      MoreArgs = list(key = key, ...),
      SIMPLIFY = FALSE
    )

  return(plots)
}


#' Dot plot
#'
#' Plot a dot plot of feature expression across identities.
#'
#' @param data The data for the plot as Seurat object or data frame.
#' @param features The name of the features to plot.
#' @param ... Additional arguments to be passed to the `geom_point` function.
#' @param cmap The colour map to use for the plot as a vector of colours. Default is `pals::kovesi.linear_blue_95_50_c20(100)`.
#' @param limits The limits for the colour scale. Optional names can be provided for the limits, which will be used as labels if
#' `binarise` is set to `TRUE`. `NA` values for a limit will infer the limit from the data. Default is c(Low = -2.5, High = 2.5).
#' @param binarise Whether to binarise the colour scale labels, taken from names of `limits`. Default is FALSE.
#' @param scale Whether to scale the expression values. Default is TRUE.
#' @param cluster Whether to cluster the identities. Default is FALSE.
#' @param min.pct The minimum percentage of cells expressing a feature to show. Default is 0.
#' @param max.size The maximum size of the points. Default is 6.
#' @param size.type The type of size scale to use, options are "radius" and "size". Default is "radius".
#' @param labs A named list of arguments to pass to the `labs` function. Default is an empty list.
#' @param theme A named list of arguments to pass to the `theme` function. Default is an empty list.
#' @param thresholds Numeric scalar or vector of the same length as `features` specifying the lower threshold value(s)
#' for the expression data. Default is NULL, which will set all thresholds to 0.
#' @param idents The name of the identity variable in the data. Default is NULL (if `data` is a Seurat object, will use the
#' default identity variable).
#' @param assay The name of the assay to search for features in if `data` is a Seurat object. Default is NULL (if `key` is
#' also not specified, will search default assay then search other assays if features not found).
#' @param layer The layer to search for features in if `data` is a Seurat object. Default is "data".
#' @param key The key of the to search for features in if `data` is a Seurat object. Default is the key associated with the
#' specified assay.
#' @param verbose Display status messages. Default is TRUE.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # Load Seurat object
#' data <- Read10X(data.dir = "path/to/data")
#' seurat_obj <- CreateSeuratObject(counts = data)
#'
#' # Plot dot plot of feature expression
#' dot.plot(data = seurat_obj, features = c("CD3E", "CD4", "CD8A"))
#'
#' @import ggplot2
#' @importFrom dplyr group_by mutate relocate ungroup
#' @importFrom pals kovesi.linear_blue_95_50_c20
#' @importFrom scales squish
#' @importFrom Seurat FetchData Idents
#'
#' @export
dot.plot <- function(data, features, ..., cmap = pals::kovesi.linear_blue_95_50_c20(100),
                     limits = c(Low = -2.5, High = 2.5), binarise = TRUE, scale = TRUE,
                     cluster = FALSE, min.pct = 0, max.size = 6, size.type = c("radius", "size"),
                     labs = list(), theme = list(), thresholds = NULL, idents = NULL,
                     assay = NULL, layer = "data", key = Seurat::Key(data[[assay]]),
                     verbose = TRUE) {
  thresholds <- if (is.null(thresholds)) rep_len(0, length(features)) else thresholds
  if (length(thresholds) != length(features)) stop("Length of thresholds must match length of features.")

  scale.size.fun <- switch(match.arg(size.type),
    radius = ggplot2::scale_radius,
    size = ggplot2::scale_size
  )

  if (!inherits(data, "Seurat")) {
    if (verbose && !is.null(assay)) message("Ignoring assay, layer and key because data is not a Seurat object.")
    if (is.null(idents)) stop("idents must be provided when data is not a Seurat object.")
    data <- data |>
      dplyr::rename("idents" = !!rlang::sym(idents)) |>
      dplyr::relocate(idents, .after = dplyr::last_col())
  } else {
    if (verbose && !is.null(assay) && !missing("key")) message("Both assay and key are provided, key will take precedence.")
    if (is.null(idents)) {
      idents <- Seurat::Idents(data)
      data <- Seurat::FetchData(object = data, vars = paste0(key, features), layer = layer)
      data[["idents"]] <- idents
    } else {
      data <- Seurat::FetchData(object = data, vars = c(paste0(key, features), idents)) |> dplyr::rename("idents" = !!rlang::sym(idents))
    }
  }

  data[["idents"]] <- factor(data[["idents"]])
  levels <- levels(data[["idents"]])

  data <- setNames(nm = levels) |>
    lapply(function(ident) {
      subset <- data[data[["idents"]] == ident, seq_len(ncol(data) - 1), drop = FALSE]
      avg.exp <- colMeans(subset)
      pct.exp <- mapply(function(x, t) mean(x > t), x = subset, t = thresholds)
      return(list(avg_exp = avg.exp, pct_exp = pct.exp))
    })

  if (cluster) {
    mat <- do.call(rbind, lapply(data, unlist)) |> scale()
    levels <- levels[hclust(dist(mat))$order]
  } else {
    levels <- rev(levels)
  }

  data <- mapply(function(data, ident) {
    as.data.frame(data) |>
      dplyr::mutate(y = ident) |>
      tibble::rownames_to_column(var = "x")
  }, data = data, ident = names(data), SIMPLIFY = FALSE)
  data <- do.call(rbind, data) |>
    dplyr::mutate(
      x = factor(sub(pattern = key, replacement = "", x = x), levels = features),
      y = forcats::fct_relevel(y, levels),
      pct_exp = ifelse(pct_exp * 100 >= min.pct, pct_exp * 100, NA)
    )

  if (scale) {
    data <- data |>
      dplyr::group_by(x) |>
      dplyr::mutate(avg_exp = scale(avg_exp)) |>
      dplyr::ungroup()
  }

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = x, y = y, colour = avg_exp, size = pct_exp)) +
    ggplot2::geom_point(...) +
    ggplot2::scale_colour_gradientn(colours = cmap, limits = limits, oob = scales::squish, breaks = if (binarise) function(x) x else ggplot2::waiver(), labels = if (binarise) names(limits) else ggplot2::waiver()) +
    scale.size.fun(range = c(0, max.size)) +
    do.call(ggplot2::labs, modifyList(list(x = "Features", y = "Identities", colour = ifelse(scale, "Scaled expression", "Expression"), size = "% Cells"), labs, keep.null = TRUE)) +
    do.call(theme_custom, modifyList(list(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)), theme, keep.null = TRUE))

  return(plot)
}


#' Inset plot legend
#'
#' This function creates a ggplot theme that moves the legend into the plot area.
#'
#' @param location The desired location of the legend within the plot area.
#'   Possible values are "topright", "topleft", "bottomright", and "bottomleft".
#'   Default is "topright".
#' @param bg The background color of the legend. Default is NA (transparent).
#'
#' @return A `ggplot` theme that moves the legend into the plot area.
#'
#' @examples
#' # Create a ggplot object
#' p <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
#'   geom_point() +
#'   labs(title = "Iris Dataset", x = "Sepal Length", y = "Sepal Width")
#'
#' # Apply the inset.legend theme to move the legend into the plot area
#' p + inset.legend(location = "topleft", bg = "white")
#'
#' @import ggplot2
#'
#' @export
inset.legend <- function(location = c("topright", "topleft", "bottomright", "bottomleft"), bg = NA) {
  location <- match.arg(location)
  anchor <- switch(location,
    topright = c(0.99, 0.99),
    topleft = c(0.01, 0.99),
    bottomright = c(0.99, 0.01),
    bottomleft = c(0.01, 0.01)
  )
  justification <- ifelse(grepl("left", location), 0, 1)
  text.position <- ifelse(grepl("left", location), "right", "left")

  theme <- ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = anchor,
    legend.justification = round(anchor),
    legend.title = ggplot2::element_text(hjust = justification),
    legend.background = ggplot2::element_rect(fill = bg),
    legend.key = ggplot2::element_blank(),
    legend.text.position = text.position
  )

  return(theme)
}


#' Calculate Centroid Density
#'
#' Calculate the centroid of each group in a dataset using 2D kernel density estimation.
#' For each group, this function finds the point with maximum density within the group's
#' data distribution. For groups with fewer than 3 points, the centroid is calculated
#' as the simple mean of the coordinates.
#'
#' @param data A data frame containing the data points.
#' @param group.var The name of the variable used for grouping the data.
#' @param x.var The name of the variable for the x-coordinate.
#' @param y.var The name of the variable for the y-coordinate.
#'
#' @return A data frame with the centroid coordinates for each group. Contains columns
#' for the grouping variable and the calculated x and y coordinates of the centroids.
#'
#' @examples
#' # Calculate centroids for clusters in a UMAP embedding
#' centroids <- calculate.centroid.density(
#'   data = umap_data,
#'   group.var = "cluster",
#'   x.var = "umap_1",
#'   y.var = "umap_2"
#' )
#'
#' @importFrom KernSmooth bkde2D
#'
#' @export
calculate.centroid.density <- function(data, group.var, x.var, y.var) {
  groups <- unique(data[[group.var]])
  res <- lapply(groups, function(g) {
    data.group <- data[data[[group.var]] == g, , drop = FALSE]
    if (nrow(data.group) < 3) {
      df <- data.frame(
        group = g,
        x = mean(data.group[[x.var]], na.rm = TRUE),
        y = mean(data.group[[y.var]], na.rm = TRUE)
      )
    } else {
      coords <- cbind(data.group[[x.var]], data.group[[y.var]])
      x.range <- range(data.group[[x.var]], na.rm = TRUE)
      y.range <- range(data.group[[y.var]], na.rm = TRUE)
      x.expand <- diff(x.range) * 0.1
      y.expand <- diff(y.range) * 0.1
      # bandwidth using rule-of-thumb
      h <- c(bw.nrd(data.group[[x.var]]), bw.nrd(data.group[[y.var]]))
      density.est <- KernSmooth::bkde2D(coords, bandwidth = h, gridsize = c(50, 50))
      max.idx <- which(density.est$fhat == max(density.est$fhat), arr.ind = TRUE)
      max.x <- density.est$x1[max.idx[1, 1]]
      max.y <- density.est$x2[max.idx[1, 2]]
      df <- data.frame(group = g, x = max.x, y = max.y)
    }
    return(df)
  })
  result <- do.call(rbind, res)
  names(result)[1] <- group.var
  return(result)
}


#' Calculate 2 dimensional point density.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param ... Additional arguments passed to `MASS::kde2d()`.
#'
#' @return Matrix of density values.
#'
#' @importFrom MASS kde2d
#'
#' @export
calculate.point.density <- function(x, y, ...) {
  d <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, d$x)
  iy <- findInterval(y, d$y)
  ii <- cbind(ix, iy)

  return(d$z[ii])
}


#' Gene Expression Program Network Plot
#'
#' Create a network plot showing correlations between gene expression programs. The function calculates
#' Spearman correlations between programs and visualizes them as a network where nodes represent programs
#' and edges represent correlations above a specified threshold.
#'
#' @param data A matrix or data frame containing gene expression program scores with programs as columns.
#' @param usage.threshold The minimum usage threshold for programs to be considered highly active in a cell. Default is 0.1.
#' @param min.correlation The minimum correlation threshold for edges in the network. Default is 0.1.
#' @param layout The layout algorithm for the network. Default is "fr" (Fruchterman-Reingold).
#' @param node.size.range A numeric vector of length 2 specifying the range of node sizes. Default is c(3, 8).
#' @param edge.width.range A numeric vector of length 2 specifying the range of edge widths. Default is c(0.5, 2).
#' @param labs A named list of arguments to pass to the `labs` function. Default is list(title = "Gene Program Co-Expression Network").
#' @param theme A named list of arguments to pass to the `theme` function. Default is an empty list.
#' @param verbose Display status messages. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A `ggplot` object of the network visualization}
#'   \item{graph}{An `igraph` object of the filtered network}
#'   \item{correlation.matrix}{The full correlation matrix between programs}
#'   \item{filtered.matrix}{The correlation matrix after applying the threshold}
#' }
#'
#' @examples
#' # Create network plot of gene expression programs
#' gep.network.plot(data = program_scores, min.correlation = 0.2)
#'
#' # Customize appearance
#' gep.network.plot(
#'   data = program_scores,
#'   usage.threshold = 0.1,
#'   min.correlation = 0.1,
#'   layout = "kk",
#'   labs = list(title = "Custom Network Title")
#' )
#'
#' @import ggplot2
#' @importFrom igraph graph_from_adjacency_matrix degree delete_vertices vcount ecount V E
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text scale_edge_width_continuous scale_edge_alpha_continuous
#' @importFrom ggrepel geom_text_repel
#' @importFrom glue glue
#'
#' @export
gep.network.plot <- function(data, usage.threshold = 0.1, min.correlation = 0.1, layout = "fr",
                             node.size.range = c(3, 8), edge.width.range = c(0.5, 2),
                             labs = list(title = "Gene Program Co-Expression Network"),
                             theme = list(), verbose = TRUE) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a matrix or data.frame")
  }

  if (ncol(data) < 2) {
    stop("data must have at least 2 columns (gene programs)")
  }

  if (usage.threshold < 0 || usage.threshold > 1) {
    stop("usage.threshold must be between 0 and 1")
  }

  if (min.correlation < 0 || min.correlation > 1) {
    stop("min.correlation must be between 0 and 1")
  }

  # Check for missing values
  if (any(is.na(data))) {
    if (verbose) message("Missing values detected in data. Using complete observations only.")
  }

  # Calculate Spearman correlation matrix between programs
  cor.matrix <- cor(data, method = "spearman", use = "complete.obs")

  # Check for any NaN or infinite values
  if (any(is.na(cor.matrix)) || any(is.infinite(cor.matrix))) {
    warning("NaN or infinite values found in correlation matrix. Check input data.")
  }

  # Set diagonal to 0 (remove self-correlations)
  diag(cor.matrix) <- 0

  # Keep only positive correlations above threshold
  filtered.matrix <- cor.matrix
  filtered.matrix[filtered.matrix < min.correlation] <- 0

  # Check if any edges remain after filtering
  if (sum(filtered.matrix > 0) == 0) {
    warning(sprintf("No positive correlations above threshold (%.2f). Consider lowering min.correlation.", min.correlation))
    return(list(plot = NULL, graph = NULL, correlation.matrix = cor.matrix, filtered.matrix = filtered.matrix))
  }

  # Create igraph object using correlation values
  g <- igraph::graph_from_adjacency_matrix(
    filtered.matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Remove isolated nodes
  isolated <- which(igraph::degree(g) == 0)
  if (length(isolated) > 0) {
    if (verbose) message(sprintf("Removing %d isolated nodes with no connections above threshold", length(isolated)))
    g <- igraph::delete_vertices(g, isolated)
  }

  # Check if graph has any vertices left
  if (igraph::vcount(g) == 0) {
    warning("No vertices remain after filtering. Consider lowering min.correlation.")
    return(list(plot = NULL, graph = NULL, correlation.matrix = cor.matrix, filtered.matrix = filtered.matrix))
  }

  if (verbose) message(sprintf("Network contains %d nodes and %d edges", igraph::vcount(g), igraph::ecount(g)))

  # Add node attributes
  igraph::V(g)$label <- igraph::V(g)$name
  igraph::V(g)$size <- colSums(data[, igraph::V(g)$name, drop = FALSE] > usage.threshold, na.rm = TRUE)

  # Calculate dynamic breaks for edge width legend based on actual correlations
  edge.weights <- igraph::E(g)$weight
  if (length(unique(edge.weights)) > 1) {
    weight.breaks <- round(seq(min(edge.weights), max(edge.weights), length.out = 4), 2)
    weight.breaks <- unique(weight.breaks)
  } else {
    weight.breaks <- round(edge.weights[1], 2)
  }

  # Create plot with refined aesthetics
  plot <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight, alpha = weight),
      color = "#34495E",
      end_cap = ggraph::circle(3, "mm"),
      start_cap = ggraph::circle(3, "mm")
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = size),
      color = "#2C3E50",
      fill = "#ECF0F1",
      shape = 21,
      stroke = 1.5
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label),
      size = 3.5,
      fontface = "bold",
      color = "#2C3E50",
      repel = TRUE,
      box.padding = 0.4,
      point.padding = 0.4,
      max.overlaps = Inf
    ) +
    ggraph::scale_edge_width_continuous(
      name = glue::glue("Correlation\n(R > {min.correlation})"),
      range = edge.width.range,
      breaks = weight.breaks,
      labels = weight.breaks,
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(color = "#34495E", alpha = 1),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        keywidth = 1.5,
        keyheight = 0.8
      )
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.4, 0.9),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(
      name = glue::glue("# Cells\n(Program Score > {usage.threshold})"),
      range = node.size.range,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(color = "#2C3E50", fill = "#ECF0F1", stroke = 1.5),
        title.position = "top",
        title.hjust = 0.5,
        keywidth = 1.2,
        keyheight = 1.2
      )
    ) +
    do.call(ggplot2::labs, labs) +
    ggplot2::theme_void() +
    do.call(ggplot2::theme, modifyList(list(
      plot.title = ggplot2::element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        color = "#2C3E50",
        margin = ggplot2::margin(b = 20)
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = ggplot2::unit(1, "cm"),
      legend.title = ggplot2::element_text(
        size = 11,
        face = "bold",
        color = "#2C3E50"
      ),
      legend.text = ggplot2::element_text(
        size = 10,
        color = "#2C3E50"
      ),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    ), theme, keep.null = TRUE))

  return(list(
    plot = plot,
    graph = g,
    correlation.matrix = cor.matrix,
    filtered.matrix = filtered.matrix
  ))
}

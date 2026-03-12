######### Diagnostics and Visualisation #####################################
#
# Tools for checking imputation validity and visualising results.

# -------------------------------------------------------------------------
# Imputation checks
# -------------------------------------------------------------------------

#' Check that all imputations fall below their censoring limits
#'
#' @param Z_orig  Original data matrix with NA for censored entries
#' @param Z_imp   A single completed data matrix (no NAs)
#' @param LOD     Censoring limits: vector (length q) or matrix (n x q)
#' @param tol     Numerical tolerance (default machine epsilon^0.5)
#' @return A list with violation summaries
check_imputation_LOD <- function(Z_orig, Z_imp, LOD,
                                 tol = .Machine$double.eps^0.5) {
  stopifnot(all(dim(Z_orig) == dim(Z_imp)))
  n <- nrow(Z_orig)
  q <- ncol(Z_orig)

  lod_mat <- if (is.vector(LOD) && length(LOD) == q) {
    matrix(LOD, nrow = n, ncol = q, byrow = TRUE)
  } else if (is.matrix(LOD) && all(dim(LOD) == c(n, q))) {
    LOD
  } else {
    stop("LOD must be a vector of length q or a matrix n x q")
  }

  miss_idx <- is.na(Z_orig)
  viol_idx <- miss_idx & (Z_imp > lod_mat + tol)

  total <- sum(viol_idx, na.rm = TRUE)

  res <- list(
    any_violation        = total > 0L,
    total_violations     = total,
    violations_by_column = colSums(viol_idx, na.rm = TRUE),
    violations_matrix    = viol_idx
  )

  if (res$any_violation) {
    warning("Found ", total, " imputations > LOD.")
  } else {
    message("All imputations are <= LOD.")
  }

  res
}

# -------------------------------------------------------------------------
# Posterior summaries and plotting helpers
# -------------------------------------------------------------------------

make_lod_matrix <- function(LOD, Z) {
  if (is.null(LOD)) {
    return(NULL)
  }

  n <- nrow(Z)
  q <- ncol(Z)

  if (is.vector(LOD) && length(LOD) == q) {
    return(matrix(LOD, nrow = n, ncol = q, byrow = TRUE))
  }

  if (is.matrix(LOD) && all(dim(LOD) == c(n, q))) {
    return(LOD)
  }

  stop("LOD must be NULL, a vector of length q, or a matrix n x q")
}

lod_reference_value <- function(lod_mat, is_imputed, j) {
  if (is.null(lod_mat)) {
    return(NA_real_)
  }

  imp_vals <- lod_mat[is_imputed[, j], j]
  imp_vals <- imp_vals[is.finite(imp_vals)]

  if (length(imp_vals) > 0L) {
    return(stats::median(imp_vals))
  }

  vals <- lod_mat[, j]
  vals <- vals[is.finite(vals)]

  if (length(vals) > 0L) {
    return(stats::median(vals))
  }

  NA_real_
}

safe_pretty_breaks <- function(xlim, n = 12) {
  br <- pretty(xlim, n = n)

  if (length(br) < 2L || any(!is.finite(br))) {
    br <- seq(xlim[1], xlim[2], length.out = n + 1L)
  }

  if (length(unique(br)) < 2L) {
    span <- max(abs(xlim), 1)
    br <- seq(xlim[1] - 0.05 * span, xlim[2] + 0.05 * span,
              length.out = n + 1L)
  }

  br
}

empty_panel <- function(bg_panel) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  rect(0, 0, 1, 1, col = bg_panel, border = NA)
}

point_alpha <- function(scores, uncertainty, alpha_range) {
  if (!(uncertainty %in% c("alpha", "both"))) {
    return(rep(alpha_range[2], length(scores)))
  }

  alpha_range[2] - scores * diff(alpha_range)
}

point_size <- function(scores, uncertainty, cex_base, cex_range) {
  if (!(uncertainty %in% c("size", "both"))) {
    return(rep(cex_base, length(scores)))
  }

  cex_range[1] + scores * diff(cex_range)
}

normalise_panel_uncertainty <- function(ps, ranges, x, y,
                                        uncertainty_stat = c("ci", "sd")) {
  uncertainty_stat <- match.arg(uncertainty_stat)

  scale_x <- diff(ranges[[x]])
  scale_y <- diff(ranges[[y]])
  if (!is.finite(scale_x) || scale_x <= 0) scale_x <- 1
  if (!is.finite(scale_y) || scale_y <= 0) scale_y <- 1

  if (uncertainty_stat == "ci") {
    ux <- ifelse(ps$is_imputed[, x], (ps$hi[, x] - ps$lo[, x]) / scale_x, 0)
    uy <- ifelse(ps$is_imputed[, y], (ps$hi[, y] - ps$lo[, y]) / scale_y, 0)
  } else {
    ux <- ifelse(ps$is_imputed[, x], ps$sd[, x] / scale_x, 0)
    uy <- ifelse(ps$is_imputed[, y], ps$sd[, y] / scale_y, 0)
  }

  scores <- sqrt(ux^2 + uy^2)
  positive <- scores[scores > 0]

  if (!length(positive)) {
    return(rep(0, length(scores)))
  }

  clip <- unname(stats::quantile(positive, probs = 0.95, na.rm = TRUE))
  if (!is.finite(clip) || clip <= 0) {
    clip <- max(positive)
  }

  pmin(scores / clip, 1)
}

representative_completed_data <- function(Y, Z, ps, burn = 0,
                                          position = c("sample", "median"),
                                          draw = NULL) {
  position <- match.arg(position)
  R <- dim(Y)[3]
  keep <- seq.int(burn + 1L, R)

  if (!length(keep)) {
    stop("No post-burn-in draws available.")
  }

  if (position == "median") {
    return(ps$median)
  }

  if (is.null(draw)) {
    draw <- keep[ceiling(length(keep) / 2)]
  }

  if (!(draw %in% keep)) {
    stop("draw must be an integer between burn + 1 and the number of draws.")
  }

  completed <- Y[, , draw]
  completed[!ps$is_imputed] <- Z[!ps$is_imputed]
  completed
}

draw_uncertainty_legend <- function(uncertainty, cex_base, cex_range,
                                    alpha_range, col_ref, bg_panel) {
  if (uncertainty == "none") {
    text(0.08, 0.86, "Uncertainty", adj = c(0, 1), font = 2, cex = 0.95)
    return(invisible(NULL))
  }

  text(0.08, 0.9, "Uncertainty", adj = c(0, 1), font = 2, cex = 0.95)

  ys <- c(0.68, 0.48, 0.28)
  xs <- rep(0.18, length(ys))
  labels <- c("Low", "Medium", "High")
  scores <- c(0.05, 0.5, 1)

  if (uncertainty == "whisker") {
    spans <- c(0.05, 0.11, 0.18)
    for (k in seq_along(scores)) {
      col_k <- grDevices::adjustcolor(col_ref, alpha.f = 0.35)
      segments(xs[k] - spans[k], ys[k], xs[k] + spans[k], ys[k],
               col = col_k, lwd = 0.8)
      segments(xs[k], ys[k] - 0.7 * spans[k], xs[k], ys[k] + 0.7 * spans[k],
               col = col_k, lwd = 0.8)
      points(xs[k], ys[k], pch = 16, cex = cex_base,
             col = grDevices::adjustcolor(col_ref, alpha.f = 0.9))
    }

  } else {
    cex_vals <- point_size(scores, uncertainty, cex_base, cex_range)
    alpha_vals <- point_alpha(scores, uncertainty, alpha_range)

    for (k in seq_along(scores)) {
      points(xs[k], ys[k], pch = 16, cex = cex_vals[k],
             col = grDevices::adjustcolor(col_ref, alpha.f = alpha_vals[k]))
    }

  }

  text(0.33, ys, labels = labels, adj = c(0, 0.5), cex = 0.8)
}

#' Compute pointwise posterior summaries for imputed values
#'
#' @param Y      Array n x q x R of posterior samples (from gibbs_mvn or
#'               gibbs_DPML output)
#' @param Z      Original data matrix (n x q, NA = censored)
#' @param alpha  Significance level for credible intervals (default 0.05)
#' @param burn   Number of initial iterations to discard (default 0)
#' @return A list with components:
#'   \item{median}{n x q matrix of pointwise posterior medians}
#'   \item{lo}{n x q matrix of lower CI bounds}
#'   \item{hi}{n x q matrix of upper CI bounds}
#'   \item{sd}{n x q matrix of posterior standard deviations}
#'   \item{is_imputed}{n x q logical matrix (TRUE = was censored)}
posterior_summary <- function(Y, Z, alpha = 0.05, burn = 0) {
  n <- nrow(Z)
  q <- ncol(Z)
  R <- dim(Y)[3]

  idx <- seq.int(burn + 1L, R)
  if (length(idx) < 2L) {
    stop("Not enough post-burn-in samples.")
  }

  med <- lo <- hi <- sd_mat <- matrix(NA_real_, n, q)
  is_imp <- is.na(Z)
  probs <- c(alpha / 2, 0.5, 1 - alpha / 2)

  for (i in seq_len(n)) {
    for (j in seq_len(q)) {
      if (is_imp[i, j]) {
        draws <- Y[i, j, idx]
        qs <- stats::quantile(draws, probs = probs, names = FALSE)
        lo[i, j] <- qs[1]
        med[i, j] <- qs[2]
        hi[i, j] <- qs[3]
        sd_mat[i, j] <- stats::sd(draws)
      } else {
        med[i, j] <- Z[i, j]
        lo[i, j] <- Z[i, j]
        hi[i, j] <- Z[i, j]
        sd_mat[i, j] <- 0
      }
    }
  }

  list(
    median     = med,
    lo         = lo,
    hi         = hi,
    sd         = sd_mat,
    is_imputed = is_imp
  )
}

#' Scatterplot matrix for multivariate LOD imputations
#'
#' The diagonal shows histograms of observed values and values below the
#' LOD. The upper triangle shows scatterplots of a single completed
#' dataset (or posterior medians), separating fully observed points from
#' the three imputed combinations. Posterior uncertainty can be encoded
#' through point alpha, point size, both, or pointwise whiskers.
#'
#' @param Y                Array n x q x R of posterior samples
#' @param Z                Original data matrix (n x q, NA = censored)
#' @param LOD              Censoring limits: NULL, vector length q, or
#'                         matrix n x q
#' @param alpha            Significance level for credible intervals
#' @param burn             Burn-in iterations to discard
#' @param varnames         Character vector of variable names
#' @param position         "sample" for one completed draw, "median" for
#'                         posterior medians
#' @param draw             Draw index to display when position = "sample"
#' @param uncertainty      Uncertainty encoding: "none", "alpha", "size",
#'                         "both", or "whisker"
#' @param uncertainty_stat Use CI width ("ci") or posterior SD ("sd") as the
#'                         uncertainty score
#' @param col_obs          Colour for fully observed points
#' @param col_imp          Optional common colour for all imputed points
#' @param col_imp_x        Colour for x-imputed / y-observed points
#' @param col_imp_y        Colour for y-imputed / x-observed points
#' @param col_imp_both     Colour for observations imputed on both axes
#' @param col_hist_obs     Histogram colour for observed values
#' @param col_hist_imp     Histogram colour for values below the LOD
#' @param col_lod          Colour for LOD reference lines
#' @param col_trend        Colour for the trend line
#' @param col_ci           Optional colour override for whiskers
#' @param pch_obs          Point character for observed values
#' @param pch_imp          Point character for imputed values
#' @param cex_obs          Point size for observed values
#' @param cex_imp          Baseline point size for imputed values
#' @param cex_range        Range used when uncertainty changes point size
#' @param alpha_range      Range used when uncertainty changes transparency
#' @param whisker_lwd      Line width for whiskers
#' @param hist_breaks      Number of pretty histogram bins
#' @param show_trend       Logical; add a least-squares line to scatterplots
#' @param show_legend      Logical; use lower-triangle panels for legends
#' @param show_varnames    Logical; print variable names on the diagonal
#' @param bg_panel         Panel background colour
#' @param main             Overall title
#' @param ...              Additional arguments passed to par()
#' @return Invisible list with posterior summary, plotted coordinates, and
#'         the draw used when position = "sample"
scatterplot_matrix_ci <- function(Y, Z, LOD = NULL, alpha = 0.05, burn = 0,
                                  varnames = colnames(Z),
                                  position = c("sample", "median"),
                                  draw = NULL,
                                  uncertainty = c("both", "alpha", "size",
                                                  "whisker", "none"),
                                  uncertainty_stat = c("ci", "sd"),
                                  col_obs = "#4B1369",
                                  col_imp = NULL,
                                  col_imp_x = "#D9A400",
                                  col_imp_y = "#F17C1E",
                                  col_imp_both = "#D62839",
                                  col_hist_obs = "#8E63A9",
                                  col_hist_imp = "#E3C34B",
                                  col_lod = "#4D4D4D",
                                  col_trend = "#2C3E99",
                                  col_ci = NULL,
                                  pch_obs = 16,
                                  pch_imp = 16,
                                  cex_obs = 0.55,
                                  cex_imp = 0.75,
                                  cex_range = c(0.75, 1.8),
                                  alpha_range = c(0.2, 0.95),
                                  whisker_lwd = 0.75,
                                  hist_breaks = 12,
                                  show_trend = TRUE,
                                  show_legend = TRUE,
                                  show_varnames = TRUE,
                                  bg_panel = "#F2F2F2",
                                  main = NULL, ...) {

  position <- match.arg(position)
  uncertainty <- match.arg(uncertainty)
  uncertainty_stat <- match.arg(uncertainty_stat)

  n <- nrow(Z)
  q <- ncol(Z)

  if (is.null(varnames)) {
    varnames <- paste0("V", seq_len(q))
  }

  if (!is.null(col_imp)) {
    col_imp_x <- col_imp_y <- col_imp_both <- col_imp
    col_hist_imp <- col_imp
  }

  ps <- posterior_summary(Y, Z, alpha = alpha, burn = burn)
  completed <- representative_completed_data(
    Y = Y, Z = Z, ps = ps, burn = burn, position = position, draw = draw
  )
  lod_mat <- make_lod_matrix(LOD, Z)

  ranges <- lapply(seq_len(q), function(j) {
    vals <- c(
      completed[, j],
      ps$lo[ps$is_imputed[, j], j],
      ps$hi[ps$is_imputed[, j], j]
    )

    if (!is.null(lod_mat)) {
      vals <- c(vals, lod_mat[, j])
    }

    rng <- range(vals, finite = TRUE, na.rm = TRUE)
    if (!all(is.finite(rng))) {
      rng <- c(0, 1)
    }

    if (diff(rng) <= 0) {
      pad <- max(abs(rng[1]), 1) * 0.05
      rng <- rng + c(-pad, pad)
    }

    pad <- 0.04 * diff(rng)
    rng + c(-pad, pad)
  })

  hist_legend_pos <- if (show_legend && q >= 3L) c(max(2L, q - 1L), 1L) else NULL
  scatter_legend_pos <- if (show_legend && q >= 2L) c(q, 1L) else NULL
  uncertainty_legend_pos <- if (show_legend && q >= 4L) c(q, 2L) else NULL

  same_panel <- function(row, col, pos) {
    !is.null(pos) && row == pos[1] && col == pos[2]
  }

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  par(
    mfrow = c(q, q),
    mar = c(1.5, 1.5, 0.9, 0.5),
    oma = c(2.5, 2.8, if (!is.null(main)) 3 else 1.5, 1),
    mgp = c(1.15, 0.25, 0),
    tcl = -0.2,
    cex.axis = 0.55,
    xaxs = "i",
    yaxs = "i",
    ...
  )

  for (row in seq_len(q)) {
    for (col in seq_len(q)) {
      if (row == col) {
        obs_vals <- completed[!ps$is_imputed[, col], col]
        imp_vals <- completed[ps$is_imputed[, col], col]
        lod_val <- lod_reference_value(lod_mat, ps$is_imputed, col)
        br <- safe_pretty_breaks(ranges[[col]], n = hist_breaks)

        h_obs <- if (length(obs_vals) > 0L) {
          hist(obs_vals, breaks = br, plot = FALSE)
        } else {
          NULL
        }

        h_imp <- if (length(imp_vals) > 0L) {
          hist(imp_vals, breaks = br, plot = FALSE)
        } else {
          NULL
        }

        ymax <- max(c(
          0,
          if (!is.null(h_obs)) h_obs$density else 0,
          if (!is.null(h_imp)) h_imp$density else 0
        ))
        if (!is.finite(ymax) || ymax <= 0) ymax <- 1

        plot.new()
        plot.window(xlim = ranges[[col]], ylim = c(0, 1.05 * ymax))
        rect(ranges[[col]][1], 0, ranges[[col]][2], 1.05 * ymax,
             col = bg_panel, border = NA)

        if (!is.null(h_imp)) {
          rect(
            xleft = h_imp$breaks[-length(h_imp$breaks)],
            ybottom = 0,
            xright = h_imp$breaks[-1],
            ytop = h_imp$density,
            col = grDevices::adjustcolor(col_hist_imp, alpha.f = 0.95),
            border = "white", lwd = 0.5
          )
        }

        if (!is.null(h_obs)) {
          rect(
            xleft = h_obs$breaks[-length(h_obs$breaks)],
            ybottom = 0,
            xright = h_obs$breaks[-1],
            ytop = h_obs$density,
            col = grDevices::adjustcolor(col_hist_obs, alpha.f = 0.95),
            border = "white", lwd = 0.5
          )
        }

        if (is.finite(lod_val)) {
          abline(v = lod_val, lty = 2, lwd = 1.1, col = col_lod)
        }

        axis(1)
        axis(2)
        if (show_varnames) {
          mtext(varnames[col], side = 3, line = 0.15, cex = 0.76, font = 2)
        }
        box(col = grDevices::adjustcolor("grey30", alpha.f = 0.55))
      } else if (row < col) {
        x <- completed[, col]
        y <- completed[, row]
        x_imp <- ps$is_imputed[, col]
        y_imp <- ps$is_imputed[, row]

        idx_obs <- !x_imp & !y_imp
        idx_x <- x_imp & !y_imp
        idx_y <- !x_imp & y_imp
        idx_both <- x_imp & y_imp

        x_lod <- lod_reference_value(lod_mat, ps$is_imputed, col)
        y_lod <- lod_reference_value(lod_mat, ps$is_imputed, row)

        scores <- normalise_panel_uncertainty(
          ps = ps, ranges = ranges, x = col, y = row,
          uncertainty_stat = uncertainty_stat
        )

        plot.new()
        plot.window(xlim = ranges[[col]], ylim = ranges[[row]])
        rect(ranges[[col]][1], ranges[[row]][1],
             ranges[[col]][2], ranges[[row]][2],
             col = bg_panel, border = NA)

        if (is.finite(x_lod)) {
          abline(v = x_lod, lty = 2, lwd = 1.1, col = col_lod)
        }
        if (is.finite(y_lod)) {
          abline(h = y_lod, lty = 2, lwd = 1.1, col = col_lod)
        }

        if (show_trend && stats::sd(x) > 0) {
          fit <- stats::lm(y ~ x)
          abline(fit, col = grDevices::adjustcolor(col_trend, alpha.f = 0.95),
                 lwd = 1.4)
        }

        if (any(idx_obs)) {
          points(x[idx_obs], y[idx_obs], pch = pch_obs, cex = cex_obs,
                 col = grDevices::adjustcolor(col_obs, alpha.f = 0.95))
        }

        draw_category <- function(idx, base_col) {
          if (!any(idx)) {
            return(invisible(NULL))
          }

          idx_num <- which(idx)
          alpha_vals <- point_alpha(scores[idx_num], uncertainty, alpha_range)
          cex_vals <- point_size(scores[idx_num], uncertainty,
                                 cex_base = cex_imp, cex_range = cex_range)

          if (uncertainty == "whisker") {
            whisker_col <- if (!is.null(col_ci)) {
              col_ci
            } else {
              grDevices::adjustcolor(base_col, alpha.f = 0.35)
            }

            idx_x_ci <- idx_num[ps$is_imputed[idx_num, col]]
            idx_y_ci <- idx_num[ps$is_imputed[idx_num, row]]

            if (length(idx_x_ci)) {
              segments(ps$lo[idx_x_ci, col], y[idx_x_ci],
                       ps$hi[idx_x_ci, col], y[idx_x_ci],
                       col = whisker_col, lwd = whisker_lwd)
            }

            if (length(idx_y_ci)) {
              segments(x[idx_y_ci], ps$lo[idx_y_ci, row],
                       x[idx_y_ci], ps$hi[idx_y_ci, row],
                       col = whisker_col, lwd = whisker_lwd)
            }
          }

          cols <- vapply(
            alpha_vals,
            function(a) grDevices::adjustcolor(base_col, alpha.f = a),
            FUN.VALUE = character(1)
          )
          points(x[idx_num], y[idx_num], pch = pch_imp, cex = cex_vals,
                 col = cols)
        }

        draw_category(idx_x, col_imp_x)
        draw_category(idx_y, col_imp_y)
        draw_category(idx_both, col_imp_both)

        axis(1)
        axis(2)

        box(col = grDevices::adjustcolor("grey30", alpha.f = 0.55))
      } else {
        empty_panel(bg_panel)

        if (same_panel(row, col, hist_legend_pos)) {
          legend(
            "center",
            legend = c("Observed values", "Values below LOD"),
            title = "Histogram legend",
            fill = c(col_hist_obs, col_hist_imp),
            border = NA,
            bty = "n",
            cex = 0.8,
            title.adj = 0
          )
        } else if (same_panel(row, col, scatter_legend_pos)) {
          legend(
            "center",
            legend = c(
              "Unobserved values due to LODs",
              "Values below LOD[x], above LOD[y]",
              "Values below LOD[y], above LOD[x]",
              "Observed values"
            ),
            title = "Scatterplot legend",
            pch = rep(pch_imp, 4),
            col = c(col_imp_both, col_imp_x, col_imp_y, col_obs),
            pt.cex = c(0.9, 0.9, 0.9, 0.8),
            bty = "n",
            cex = 0.72,
            title.adj = 0
          )
        } else if (same_panel(row, col, uncertainty_legend_pos)) {
          draw_uncertainty_legend(
            uncertainty = uncertainty,
            cex_base = cex_imp,
            cex_range = cex_range,
            alpha_range = alpha_range,
            col_ref = col_imp_both,
            bg_panel = bg_panel
          )
        }
      }
    }
  }

  if (!is.null(main)) {
    mtext(main, outer = TRUE, side = 3, line = 1.1, cex = 1.1, font = 2)
  }

  if (!show_legend || q < 4L) {
    note <- switch(
      uncertainty,
      none = "Observed and imputed values shown separately.",
      alpha = "Uncertainty encoding: alpha.",
      size = "Uncertainty encoding: size.",
      both = "Uncertainty encoding: alpha and size.",
      whisker = paste0("Uncertainty encoding: ", round(100 * (1 - alpha)),
                       "% credible interval whiskers.")
    )

    mtext(note, outer = TRUE, side = 1, line = 0.6, cex = 0.72, col = "grey35")
  }

  invisible(list(
    summary = ps,
    plotted = completed,
    draw = if (position == "sample") {
      if (is.null(draw)) seq.int(burn + 1L, dim(Y)[3])[ceiling((dim(Y)[3] - burn) / 2)]
      else draw
    } else {
      NA_integer_
    }
  ))
}

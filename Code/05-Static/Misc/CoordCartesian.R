UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 3,
  panel_counter = 1,
  layout = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    self$layout <- layout # store for later
    layout
  },
  
  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    train_cartesian <- function(scale, limits, name, given_range = c(NA, NA)) {
      if (anyNA(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion, coord_limits = limits)
        isna <- is.na(given_range)
        given_range[isna] <- range[isna]
      }
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, given_range),
        sec = ggplot2:::view_scale_secondary(scale, limits, given_range),
        arrange = scale$axis_order(),
        range = given_range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    this_layout <- self$layout[ self$panel_counter,, drop = FALSE ]
    self$panel_counter <- 
      if (self$panel_counter < self$num_of_panels) {
        self$panel_counter + 1
      } else 1
    
    # determine merge column names by removing all "standard" names
    layout_names <- setdiff(names(this_layout),
                            c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
    limits_names <- setdiff(names(self$panel_limits),
                            c("xmin", "xmax", "ymin", "ymax"))
    
    limit_extras <- setdiff(limits_names, layout_names)
    if (length(limit_extras) > 0) {
      stop("facet names in 'panel_limits' not found in 'layout': ",
           paste(sQuote(limit_extras), collapse = ","))
    } else if (length(limits_names) == 0 && NROW(self$panel_limits) == 1) {
      # no panels in 'panel_limits'
      this_panel_limits <- cbind(this_layout, self$panel_limits)
    } else {
      this_panel_limits <- merge(this_layout, self$panel_limits, all.x = TRUE, by = limits_names)
    }
    
    if (isTRUE(NROW(this_panel_limits) > 1)) {
      stop("multiple matches for current panel in 'panel_limits'")
    }
    
    # add missing min/max columns, default to "no override" (NA)
    this_panel_limits[, setdiff(c("xmin", "xmax", "ymin", "ymax"),
                                names(this_panel_limits)) ] <- NA
    
    c(train_cartesian(scale_x, self$limits$x, "x",
                      unlist(this_panel_limits[, c("xmin", "xmax"), drop = TRUE])),
      train_cartesian(scale_y, self$limits$y, "y",
                      unlist(this_panel_limits[, c("ymin", "ymax"), drop = TRUE])))
  }
)

coord_cartesian_panels <- function(panel_limits, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords,
                   panel_limits = panel_limits,
                   expand = expand, default = default, clip = clip)
}
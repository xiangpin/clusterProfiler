
#' @method plot interpretation
#' @title plot
#' @param x An `interpretation` object.
#' @param layout Graph layout, default is "nicely".
#' @param ... Additional arguments passed to `ggtangle::ggtangle`.
#' @export
plot.interpretation <- function(x, layout = "nicely", ...) {
    
    if (is.null(x$network)) {
        message("No refined network available to plot.")
        return(invisible(NULL))
    }
    
    if (!requireNamespace("ggtangle", quietly = TRUE)) {
        message("Package 'ggtangle' is required for this plot. Please install it.")
        return(invisible(NULL))
    }
    
    g <- x$network
    
    # Use ggtangle to visualize the network
    # Note: ggtangle extends ggplot2
    p <- ggtangle::ggtangle(g, layout = layout, ...)
    
    # Add edges
    # We try to map interaction type to color if available
    if ("interaction" %in% igraph::edge_attr_names(g)) {
        p <- p + ggtangle::geom_edge(ggplot2::aes(color = interaction)) +
             ggplot2::scale_color_manual(values = c(
                 activation = "green3", 
                 inhibition = "red3", 
                 repression = "red3",
                 binding = "blue",
                 interaction = "grey50"
             ), na.value = "grey50")
    } else {
        p <- p + ggtangle::geom_edge(color = "grey50")
    }
    
    # Add nodes and labels
    p <- p + ggtangle::geom_node_point(size = 5, color = "lightblue") + 
         ggtangle::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
         ggplot2::theme_void()
    
    main_title <- "Refined Regulatory Network"
    if (!is.null(x$cluster)) {
        main_title <- paste0(main_title, " (", x$cluster, ")")
    }
    
    p <- p + ggplot2::ggtitle(main_title)
    
    return(p)
}

#' @method plot interpretation_list
#' @export
plot.interpretation_list <- function(x, ...) {
    plots <- lapply(x, function(res) plot(res, ...))
    # Filter nulls
    plots <- plots[!sapply(plots, is.null)]
    
    if (length(plots) == 0) return(invisible(NULL))
    
    # Print plots
    for (p in plots) {
        print(p)
        if (length(plots) > 1 && interactive()) {
            readline(prompt="Press [enter] to see the next plot...")
        }
    }
    invisible(plots)
}

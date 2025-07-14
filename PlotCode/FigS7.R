make.grid.plots = function(min_range, max_range, by, n = NA) {

  library(cowplot)
  library(grid)

  # Create combinations of upper and lower where upper >= lower
  thresh_mat = expand.grid(lower = seq(min_range[1], min_range[2], by = by),
                           upper = seq(max_range[1], max_range[2], by = by))
  thresh_mat = thresh_mat[thresh_mat$upper >= thresh_mat$lower, ]

  # Add the sigma column
  thresh_mat$sigma = ifelse(thresh_mat$upper == thresh_mat$lower, 0,
                            (thresh_mat$upper - thresh_mat$lower) / 3)

  if (is.na(n)) {n = length(thresh_mat[, 1])}

  ### Find successful runs #####################################################
  files = list.files("../../data/Reviewer_Grid_Runs/MIN_GLU")

  runs = as.numeric(regmatches(files, regexpr("[0-9]+(?=\\.RDS$)",
                                              files, perl = TRUE)))

  ### Organize Trackers ########################################################

  trackers_list = vector("list", n)

  for (i in runs) {
    track = readRDS(paste0(
      "../../data/Reviewer_Grid_Runs/POST_TRACK/Post_track_", i, ".RDS"))
    trackers_list[[i]] = track
  }

  ### Fancy Tracker Trajectory Plot ##############################################

  line_plots_list = vector("list", n)

  for (i in 1:n) {
    if (i %in% runs) {
      line_plots_list[[i]] = ggplot(trackers_list[[i]][2001:2200, ]) +
        geom_line(aes(x = Time - 2000, y = Out_tht_sig), color = "#F8766D") +
        scale_y_continuous(n.breaks = 3, limits = c(0, 0.65)) +
        scale_x_continuous(n.breaks = 3) +
        labs(
          # title = paste0("[", c(thresh_mat[i, 1]), ", ",
          #               c(thresh_mat[i, 2]), "]"),
          x = "", # "Time",
          y = ""
        ) +
        theme(
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"),
          axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 11)
        )
    } else {
      line_plots_list[[i]] = ggplot() + theme_void()
    }
  }

  # Create a blank grob to represent empty spaces
  blank_grob = grid.rect(gp = gpar(col = "white"))

  # Prepare a 9x9 matrix of grobs
  plot_matrix = matrix(list(blank_grob),
                       nrow = round((max_range[2] - max_range[1])/by) + 1,
                       ncol = round((min_range[2] - min_range[1])/by) + 1)

  # Place plots in their specified coordinates
  for (i in 1:n) {
    x = round((thresh_mat[i, 1] - min_range[1]) * (1/by) + 1)
    y = round((max_range[2] - max_range[1])*(1/by) -
                ((thresh_mat[i, 2] - max_range[1]) * (1/by)) + 1)
    plot_matrix[[y, x]] = line_plots_list[[i]]
  }

  # Arrange the plots
  master_plot = grid.arrange(
    grobs = plot_matrix,
    layout_matrix = matrix(1:length(plot_matrix),
                           nrow = dim(plot_matrix)[1],
                           ncol = dim(plot_matrix)[2])
  )

  # Put plots on an underlying axis
  underlying_plot = ggplot() +
    xlim(min_range[1], min_range[2]) +
    ylim(max_range[1], max_range[2]) +
    labs(x = "Minimum Threshold", y = "Maximum Threshold") +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white"),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 22)
    )

  combined = ggdraw(underlying_plot) +
    draw_plot(master_plot,
              x = 0.07,     # x position (0-1 scale)
              y = 0.06,     # y position (0-1 scale)
              width = 0.9, # width of the overlay plot (fraction of original plot)
              height = 0.92 # height of the overlay plot (fraction of original plot)
    )

  ggsave("Figures/FigS7.tiff", combined,
         width = 15, height = 14, units = "in", dpi = 300)

  ggsave("Figures/FigS7.pdf", combined,
         width = 15, height = 14, units = "in", dpi = 300)
}

make.grid.plots(c(-0.6, 0.6), c(2.5, 3.5), 0.1)

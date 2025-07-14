data = data.frame(
  Time = 1:100,
  K_m = c(rep(8, 20), rep(300, 5), rep(8, 75)),
  G_m = rep(19.5, 100),
  Thresholds = rnorm(1000, -1, 0.03),
  Glutamate = c(rep(30, 30), rep(35, 40), rep(30, 30))
)

k_data = readRDS(paste0(plotdata_folder, "K_Synch_Data_norm.RDS"))


K_plot = ggplot() +
  geom_segment(aes(x = 0, xend = 27, y = 8, yend = 8),
               color = "gold", size = 0.6) +
  geom_segment(aes(x = 27, xend = 30, y = 30, yend = 30),
               color = "gold", size = 0.6) +
  geom_segment(aes(x = 30, xend = 100, y = 8, yend = 8),
               color = "gold", size = 0.6) +
  scale_y_continuous(breaks = c(0, 8, 30), labels = c("0", "8", "300"),
                     limits = c(0, 35)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  labs(title = NULL,
       x = "Time",
       y = expression("K"["m"]^"+"),
       linetype = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
#    text = element_text(size = 12)
  )

# ggsave("Plot_Code/Final_Plots/SchematicPlots/Increased_K.tiff", K_plot,
#        width = 2, height = 1, units = "in", dpi = 300)


G_red_plot = ggplot() +
  geom_segment(aes(x = 0, xend = 30, y = 30, yend = 30), color = "purple",
               size = 0.6) +
  geom_segment(aes(x = 30, xend = 100, y = 20, yend = 20), color = "purple",
               size = 0.6) +
  scale_y_continuous(breaks = c(0, 20, 30),
                     labels = c("0", "20", "30"),
                     limits = c(0, 35)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  labs(title = NULL,
       x = "Time",
       y = expression("G"["m"]),
       linetype = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    #    text = element_text(size = 12)
  )

# ggsave("Plot_Code/Final_Plots/SchematicPlots/Reduced_G.tiff", G_red_plot,
#        width = 2, height = 1, units = "in", dpi = 300)


Thresholds_plot = ggplot(data) +
  geom_density(aes(x = Thresholds), color = "black", fill = "slategray", alpha = 0.7) +
  labs(title = NULL,
       x = "Thresholds",
       y = "density",
       linetype = NULL) +
  scale_x_continuous(limits = c(-2, 4), breaks = c(-1, 0, 1, 3),
                     expand = c(0.01, 0.01)) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    #    text = element_text(size = 12)
  )

# ggsave("Plot_Code/Final_Plots/SchematicPlots/Threhsolds.tiff", Thresholds_plot,
#        width = 1.8, height = 1, units = "in", dpi = 300)


G_inc_plot = ggplot() +
  geom_segment(aes(x = 0, xend = 30, y = 30, yend = 30),
               color = "purple", size = 0.6) +
  geom_segment(aes(x = 30, xend = 70, y = 42, yend = 42),
               color = "purple", size = 0.6) +
  geom_segment(aes(x = 70, xend = 100, y = 30, yend = 30),
               color = "purple", size = 0.6) +
  scale_y_continuous(breaks = c(0, 30, 42), labels = c("0", "30", "35"),
                     limits = c(0, 50)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  labs(title = NULL,
       x = "Time",
       y = expression("G"["m"]),
       linetype = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    #    text = element_text(size = 12)
  )

# ggsave("Plot_Code/Final_Plots/SchematicPlots/Increased_G.tiff", G_inc_plot,
#        width = 2, height = 1, units = "in", dpi = 300)



# setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))
# plot_data = read.csv("Plot_Code/Plot_Data/Synchronization_Data_1.csv", header = TRUE)
#
# g_osc = ggplot() +
#   geom_line(aes(x = plot_data$time[151:189] - 1, y = plot_data$g_traj[151:189]),
#             color = "purple", size = 0.6) +
#   geom_line(aes(x = plot_data$time[111:150], y = plot_data$g_traj[111:150]),
#             color = "purple", size = 0.6) +
#   geom_segment(aes(x = 150, xend = 150, y = plot_data$g_traj[150],
#                    yend = plot_data$g_traj[151]),
#                lty = "12", color = "purple", size = 0.6) +
#   scale_y_continuous(breaks = c(min(plot_data$g_traj), max(plot_data$g_traj)),
#                      labels = c("25", "35"),
#                      limits = c(min(plot_data$g_traj) - 0.05,
#                                 max(plot_data$g_traj) + 0.05)) +
#   scale_x_continuous(expand = c(0.01, 0.01), limits = c(111, 188)) +
#   labs(title = NULL,
#        x = NULL,
#        y = expression("G"["m"]),
#        linetype = NULL) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "white"),
#     plot.background = element_rect(fill = "white", color = "white"),
#     axis.line = element_line(color = "black"),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#     #    text = element_text(size = 12)
#   )

k_osc = ggplot() +
  geom_line(aes(y = k_data$k_traj[461:601], x = 1:141),
            color = "gold", size = 0.6) +
  # geom_line(aes(x = plot_data$time[151:189] - 1, y = plot_data$k_traj[151:189]),
  #           color = "gold", size = 0.6) +
  # geom_line(aes(x = plot_data$time[111:150], y = plot_data$k_traj[111:150]),
  #           color = "gold", size = 0.6) +
  # geom_segment(aes(x = 150, xend = 150, y = plot_data$k_traj[150],
  #                  yend = plot_data$k_traj[151]),
  #              lty = "12", color = "gold", size = 0.6) +
  # scale_y_continuous(breaks = c(min(plot_data$k_traj), max(plot_data$k_traj)),
  #                    labels = c("7.7", "8.3"),
  #                    limits = c(min(plot_data$k_traj) - 0.05,
  #                               max(plot_data$k_traj) + 0.05)) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(1, 141)) +
  labs(title = NULL,
       x = "Time",
       y = expression("K"["m"]^"+"),
       linetype = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    #    text = element_text(size = 12)
  )

# g_osc / k_osc


# ggsave("Plot_Code/Final_Plots/SchematicPlots/Oscillation.tiff", g_osc / k_osc,
#        width = 2.3, height = 2, units = "in", dpi = 300)
#
#
# K_plot / G_red_plot / Thresholds_plot / G_inc_plot / g_osc / k_osc

ggsave("PlotCode/Plot_Data/SchematicComponent5.tiff",
       K_plot / G_red_plot / Thresholds_plot / G_inc_plot / k_osc,
       width = 2.3, height = 5, units = "in", dpi = 300)

ggsave("PlotCode/Plot_Data/SchematicComponent6.tiff",
       K_plot / G_red_plot / Thresholds_plot / G_inc_plot / k_osc,
       width = 2.3, height = 6, units = "in", dpi = 300)

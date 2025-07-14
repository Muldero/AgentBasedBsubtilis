setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

which_cell = 6

pot_data = readRDS(
  "Output/Other/Potassium_Trajectories1.RDS")[which_cell, ] %>% t %>% c
glu_data = readRDS(
  "Output/Other/Glutamate_Trajectories_1.RDS")[which_cell, ] %>% t %>% c
mem_pot_data = readRDS("Output/Other/Membrane_Potential_Trajectories_1.RDS"
)[which_cell, ] %>% t %>% c

dif = 9

mn = 800 - dif
mx = 1101 - dif

data = data.frame(
  Time = 1:(mx - mn + 1),
  Pot = pot_data[mn:mx],
  Glu = glu_data[mn:mx],
  Mem = mem_pot_data[mn:mx]
)

glu_plot = data %>% ggplot() +
  geom_line(aes(x = Time, y = Glu), color = "purple") +
  labs(title = bquote(bold(A)),
       x = element_blank(), # "Time",
       y = "Internal Glutamate",
       linetype = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

pot_plot = data %>% ggplot() +
  geom_line(aes(x = Time, y = Pot), color = "gold") +
  labs(title = bquote(bold(B)),
       x = element_blank(), # "Time",
       y = "Internal Potassium",
       linetype = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

mem_pot_plot = data %>% ggplot() +
  geom_line(aes(x = Time, y = Mem), color = "cyan") +
  scale_y_reverse() +
  labs(title = bquote(bold(C)),
       x = "Time",
       y = "Membrane Potential",
       linetype = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

cell_traces_plot = glu_plot / pot_plot / mem_pot_plot

cell_traces_plot

ggsave(paste0(return_folder, "FigS1.tiff"),
       cell_traces_plot, width = 5.2, height = 5.5, units = "in", dpi = 300)

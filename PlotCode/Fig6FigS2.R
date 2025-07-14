library(patchwork)
library(cowplot)

data_path = "../../data/Reviewer_Grid_Runs/"

### Create Phase Plot

min_range = c(-0.6, 0.6)
max_range = c(2.5, 3.5)
by = 0.1

# Create combinations of upper and lower where upper >= lower
thresh_mat = expand.grid(lower = seq(min_range[1], min_range[2], by = by),
                         upper = seq(max_range[1], max_range[2], by = by))
thresh_mat = thresh_mat[thresh_mat$upper >= thresh_mat$lower, ]

# Add the sigma column
thresh_mat$sigma = ifelse(thresh_mat$upper == thresh_mat$lower, 0,
                          (thresh_mat$upper - thresh_mat$lower) / 3)

thresh_mat = rbind(thresh_mat, c(-1, -0.1, 0.1))

n = length(thresh_mat[, 1])

trackers_list = vector("list", n - 1)

# Find successful runs
files = list.files(paste0(data_path, "MIN_GLU"))
runs = as.numeric(regmatches(files, regexpr("[0-9]+(?=\\.RDS$)",
                                            files, perl = TRUE)))

for (i in runs) {
  if (i != 144) {
    track = readRDS(paste0(data_path, "POST_TRACK/Post_track_", i, ".RDS"))
    trackers_list[[i]] = track
  }
}

################################################################################
### Make Heatmap Matrix ########################################################
################################################################################

get.range = function(vector) {
  max_peaks = get.only.max.peaks(vector, midline = median(vector))
  min_peaks = get.only.max.peaks(vector, midline = median(vector), top = FALSE)

  mean(vector[max_peaks]) - mean(vector[min_peaks])
}

heatmap_df = data.frame(x = round(thresh_mat$lower, 2),
                        y = round(thresh_mat$upper, 2),
                        value = 0,
                        value2 = 0)

heatmap_df = heatmap_df[1:(n - 1), ]

for (i in 1:(n - 1)) {
  if (i %in% runs) {
    vec = trackers_list[[i]][2401:3000, ]$Out_tht_sig
    k = trackers_list[[i]][2401:3000, ]$Out_K_e
    val = max(vec) - min(vec) # get.range(vec)
    if (val < 0.2 & mean(k) < 50) {
      heatmap_df[i, 3] = 1 #+ runif(1, -0.1, 0.1) # no sig
    } else if (mean(k) > 50) {
      heatmap_df[i, 3] = 3 #+ runif(1, -0.1, 0.1) # collapse
    } else {
      heatmap_df[i, 3] = 2 #+ runif(1, -0.1, 0.1) # stable oscillations
    }
  } else {
    heatmap_df[i, 3] = 3
  }
}

# heatmap_df$value[c(17, 72, 115)] = 5

phase_df = data.frame(
  x = c(-0.65, -0.45, -0.45, -0.15, -0.15, 0.15, 0.15, 0.45, 0.45, # 0.65,
        -0.45, -0.45, -0.25, -0.25, -0.15, -0.15, -0.05, -0.05, 0.15, 0.15,
        0.25, 0.25, 0.35, 0.35, 0.45, 0.45, 0.65),

  y = c(2.85, 2.85, 2.75, 2.75, 2.65, 2.65, 2.55, 2.55, 2.45, # 2.45,
        3.55, 3.35, 3.35, 3.25, 3.25, 3.15, 3.15, 3.05, 3.05, 2.95, 2.95,
        2.85, 2.85, 2.75, 2.75, 2.65, 2.65),

  linenum = c(rep(1, 9), rep(2, 17))
)

pp2 = ggplot() +
  geom_tile(data = heatmap_df, aes(x = x, y = y, fill = factor(value))) +
  geom_line(data = phase_df, aes(x = x, y = y, group = linenum), size = 0.6) +
  # Annotations
  annotate("text", x = -0.42, y = 2.55,
           label = "Minimal Signaling",
           color = "black", size = 4, hjust = 0) +
  annotate("text", x = -0.25, y = 2.88,
           label = "Stable Oscillations",
           color = "black", size = 4, hjust = 0) +
  annotate("text", x = 00, y = 3.3,
           label = "Biofilm Collapse",
           color = "black", size = 4, hjust = 0) +
  scale_fill_manual(
    values = c("1" = "lightgray", "2" = "slategray", "3" = "darkgray",
               "5" = "white"),
    name = "Value"
  ) +
  scale_y_continuous(breaks = c(2.5, 2.75, 3, 3.25, 3.5)) +
  scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  labs(x = "Minimum Threshold",
       y = "Maximum Threshold",
       title = bquote(bold(G)),
       fill = NULL) +
  # coord_fixed() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

pp2

################################################################################
### Make Trace Plots ###########################################################
################################################################################

# no signaling, minimal signaling (w/ inset), stable, and uncoordinated

# 000 - no signaling - row 144
# minimal signaling - -0.4, 2.5 # row 3
# stable - 031 # row 72
# uncoordinated 0.3, 2.9 # row 62

### No signaling ###############################################################
# Get data, clean, make trace and glu plot
#
# vals_nosig = data.frame(
#   "ID" = 1:length(atr[, 1]),
#   "Value" = c(readRDS("Plot_Code/Plot_Data/Efficacy_NonVals.RDS"))[[1]])
#
# traj_nosig = readRDS(
#   "Final_Pipeline/Final_Out/POST_TRACK/Post_track_Thresh[0,0]_SD0.RDS")
#
# Glu_nosig = ggplot(vals_nosig %>% filter(Value > 0)) +
#   geom_density(aes(x = Value), color = "slategray", fill = "purple",
#                alpha = 0.5) +
#   annotate("text", x = 0,
#            y = max(density((vals_nosig %>% filter(Value > 0))$Value)$y)*0.9,
#            label = "Threshold = 0",
#            color = "black", size = 3.5, hjust = 0) +
#   annotate("text", x = 0,
#            y = max(density((vals_nosig %>% filter(Value > 0))$Value)$y)*0.6,
#            label = paste0(sum(vals_nosig$Value > 0 &
#                                 vals_nosig$Value < 0.001), " starved cells"),
#            color = "black", size = 2.7, hjust = 0) +
#   labs(
#     x = "", # "Time-Averaged Internal Glutamate",
#     y = "", # "Cell Frequency",
#     title = bquote(bold(A)) # Threshold = 0"
#   ) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "white"),
#     plot.background = element_rect(fill = "white", color = "white"),
#     axis.line = element_line(color = "black"),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 11)
#   )
#
# Trac_nosig = ggplot(traj_nosig[2501:2700, ]) +
#   geom_line(aes(x = Time - 2500, y = Tot_fv_sig), color = "#F8766D") +
#   scale_y_continuous(n.breaks = 3, limits = c(0, 0.65)) +
#   labs(
#     x = "", # "Time",
#     y = "", # "Signaling Fraction",
#     title = bquote(bold(B))
#   ) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "white"),
#     plot.background = element_rect(fill = "white", color = "white"),
#     axis.line = element_line(color = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 11)
#   )

### Min signaling ##############################################################
# Get data, clean, make trace and glu plot
# -0.4, 2.5 row 3

mg = readRDS(paste0(data_path, "MIN_GLU/Minimum_Glutamate_3.RDS"))
colnames(mg) = mg[1, ]
mg = mg[700:length(mg[, 1]), ]
mean_glu = colMeans(mg)
vals = rep(0, length(atr[, 1]))
vals[as.numeric(colnames(mg))] = mean_glu

vals_minsig = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = vals)

traj_minsig = readRDS(paste0(data_path, "POST_TRACK/Post_track_3.RDS"))

Glu_minsig = ggplot(vals_minsig %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text", x = 0.2,
           y = max(density((vals_minsig %>% filter(Value > 0))$Value)$y)*0.85,
           label = paste0(sum(vals_minsig$Value > 0 &
                                vals_minsig$Value < 0.001),
                          " starved cells"),
           color = "black", size = 2.7, hjust = 0) +
  xlim(0, 3.2) +
  labs(
    x = "", # "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold('D'))# ~ "  Minimal Signaling") # Threshold = 0"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

Trac_minsig = ggplot(traj_minsig[2501:2700, ]) +
  geom_line(aes(x = Time - 2500, y = Out_tht_sig), color = "#F8766D") +
  annotate("text",
           x = 0,
           y = 0.52,
           label = "Threshold [-0.4, 2.5]",
           color = "black", size = 3.5, hjust = 0) +
  scale_y_continuous(n.breaks = 4, limits = c(0, 0.55)) +
  labs(
    x = "", # "Time",
    y = "Signaling Fraction",
    title = bquote(bold(A))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

# scale_factor_min = 0.55/max(traj_minsig$Out_K_e)
#
# Trac_minsig = ggplot(traj_minsig[2501:2700, ]) +
#   geom_line(aes(x = Time - 2500, y = Out_tht_sig),
#             color = "#F8766D") +
#   geom_line(aes(x = Time - 2500, y = Out_K_e*scale_factor_min),
#             color = "gold", lty = "solid") +
#   scale_y_continuous(
#     n.breaks = 4,
#     limits = c(0, 0.55),
#     name = "Signaling Fraction",
#     sec.axis = sec_axis(trans = ~./scale_factor_min,
#                         name = "")
#   ) +
#   annotate("text",
#            x = 0,
#            y = 0.52,
#            label = "Threshold [-0.4, 2.5]",
#            color = "black", size = 3.5, hjust = 0) +
#   # scale_y_continuous(n.breaks = 4, limits = c(0, 0.55)) +
#   labs(
#     x = "", # "Time",
#     y = "", # "Signaling Fraction"#,
#     title = bquote(bold(A))
#   ) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "white"),
#     plot.background = element_rect(fill = "white", color = "white"),
#     axis.line.y.right = element_line(color = "gold"),
#     axis.title.y.right = element_text(color = "gold"),
#     axis.text.y.right = element_text(color = "gold"),
#     axis.line.y.left = element_line(color = "#F8766D"),
#     axis.title.y.left = element_text(color = "#F8766D"),
#     axis.text.y.left = element_text(color = "#F8766D"),
#     axis.line.x = element_line(color = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 11)
#   )

### Stable signaling ###########################################################
# Get data, clean, make trace and glu plot

mg = readRDS(paste0(data_path, "MIN_GLU/Minimum_Glutamate_72.RDS"))
colnames(mg) = mg[1, ]
mg = mg[700:length(mg[, 1]), ]
mean_glu = colMeans(mg)
vals = rep(0, length(atr[, 1]))
vals[as.numeric(colnames(mg))] = mean_glu

vals_stabsig = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = vals)

traj_stabsig = readRDS(paste0(data_path, "POST_TRACK/Post_track_72.RDS"))


Glu_stabsig = ggplot(vals_stabsig %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text",
           x = 0.2,
           y = max(density((vals_stabsig %>% filter(Value > 0))$Value)$y)*0.85,
           label = paste0(sum(vals_stabsig$Value > 0 &
                                vals_stabsig$Value < 0.001),
                          " starved cells"),
           color = "black", size = 2.7, hjust = 0) +
  xlim(0, 3.2) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "", # "Cell Frequency",
    title = bquote(bold('E')) # ~ "  Stable Signaling") # Threshold = 0"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

Trac_stabsig = ggplot(traj_stabsig[2501:2700, ]) +
  geom_line(aes(x = Time - 2500, y = Out_tht_sig), color = "#F8766D") +
  annotate("text", x = 0,
           y = 0.52,
           label = "[0, 3]",
           color = "black", size = 3.5, hjust = 0) +
  scale_y_continuous(n.breaks = 4, limits = c(0, 0.55)) +
  labs(
    x = "Time",
    y = "", # "Signaling Fraction",
    title = bquote(bold(B))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

# scale_factor_stab = 0.55/max(traj_stabsig$Out_K_e)
#
# Trac_stabsig = ggplot(traj_stabsig[2501:2700, ]) +
#   geom_line(aes(x = Time - 2500, y = Out_tht_sig),
#             color = "#F8766D") +
#   geom_line(aes(x = Time - 2500, y = Out_K_e*scale_factor_stab),
#             color = "gold", lty = "solid") +
#   scale_y_continuous(
#     n.breaks = 4,
#     limits = c(0, 0.55),
#     name = "",
#     sec.axis = sec_axis(trans = ~./scale_factor_stab,
#                         name = "")
#   ) +
#   annotate("text",
#            x = 0,
#            y = 0.52,
#            label = "[0.3, 2.9]",
#            color = "black", size = 3.5, hjust = 0) +
#   labs(
#     x = "", # "Time",
#     y = "", # "Signaling Fraction"#,
#     title = bquote(bold(B))
#   ) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "white"),
#     plot.background = element_rect(fill = "white", color = "white"),
#     axis.line.y.right = element_line(color = "gold"),
#     axis.title.y.right = element_text(color = "gold"),
#     axis.text.y.right = element_text(color = "gold"),
#     axis.line.y.left = element_line(color = "#F8766D"),
#     axis.title.y.left = element_text(color = "#F8766D"),
#     axis.text.y.left = element_text(color = "#F8766D"),
#     axis.line.x = element_line(color = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 11)
#   )


### Uncoordinated signaling ####################################################
# Get data, clean, make trace and glu plot
# 0.3, 2.9 row 62

mg = readRDS(paste0(data_path, "MIN_GLU/Minimum_Glutamate_62.RDS"))
colnames(mg) = mg[1, ]
mg = mg[700:length(mg[, 1]), ]
mean_glu = colMeans(mg)
vals = rep(0, length(atr[, 1]))
vals[as.numeric(colnames(mg))] = mean_glu

vals_uncsig = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = vals)

traj_uncsig = readRDS(paste0(data_path, "POST_TRACK/Post_track_62.RDS"))


Glu_uncsig = ggplot(vals_uncsig %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text", x = 0.2,
           y = 3.3*10^30*0.85,
             # max(density((vals_uncsig %>% filter(Value > 0))$Value)$y)*0.85,
           label = paste0(sum(vals_uncsig$Value > 0 &
                                vals_uncsig$Value < 0.001),
                          " starved cells"),
           color = "black", size = 2.7, hjust = 0) +
  xlim(0, 3.2) +
  labs(
    x = "", # "Time-Averaged Internal Glutamate",
    y = "", # "Cell Frequency",
    title = bquote(bold('F')) # ~ "  Uncoordinated Signaling") # Threshold = 0"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

Trac_uncsig = ggplot(traj_uncsig[2501:2700, ]) +
  geom_line(aes(x = Time - 2500, y = Out_tht_sig), color = "#F8766D") +
  annotate("text", x = 0,
           y = 0.52,
           label = "[0.3, 2.9]",
           color = "black", size = 3.5, hjust = 0) +
  scale_y_continuous(n.breaks = 4, limits = c(0, 0.55)) +
  labs(
    x = "Time",
    y = "", # "Signaling Fraction",
    title = bquote(bold(C))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

# scale_factor_unc = 0.55/max(traj_uncsig$Out_K_e)
#
# Trac_uncsig = ggplot(traj_uncsig[2501:2700, ]) +
#   geom_line(aes(x = Time - 2500, y = Out_tht_sig),
#             color = "#F8766D") +
#   geom_line(aes(x = Time - 2500, y = Out_K_e*scale_factor_unc),
#             color = "gold", lty = "solid") +
#   scale_y_continuous(
#     n.breaks = 4,
#     limits = c(0, 0.55),
#     name = "",
#     sec.axis = sec_axis(trans = ~./scale_factor_unc,
#                         name = "Potassium (mM)")
#   ) +
#   annotate("text",
#            x = 0,
#            y = 0.52,
#            label = "[0.3, 2.9]",
#            color = "black", size = 3.5, hjust = 0) +
#   # scale_y_continuous(n.breaks = 4, limits = c(0, 0.55)) +
#   labs(
#     x = "", # "Time",
#     y = "", # "Signaling Fraction"#,
#     title = bquote(bold('C'))
#   ) +
#   theme(
#     panel.background = element_rect(fill = "white", color = "white"),
#     plot.background = element_rect(fill = "white", color = "white"),
#     axis.line.y.right = element_line(color = "gold"),
#     axis.title.y.right = element_text(color = "gold"),
#     axis.text.y.right = element_text(color = "gold"),
#     axis.line.y.left = element_line(color = "#F8766D"),
#     axis.title.y.left = element_text(color = "#F8766D"),
#     axis.text.y.left = element_text(color = "#F8766D"),
#     axis.line.x = element_line(color = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 11)
#   )
#
# Trac_uncsig

################################################################################
### Make atr plots #############################################################
################################################################################

mg = readRDS(paste0(data_path, "MIN_GLU/Minimum_Glutamate_144.RDS"))
colnames(mg) = mg[1, ]
mg = mg[700:length(mg[, 1]), ]
mean_glu = colMeans(mg)
vals = rep(0, length(atr[, 1]))
vals[as.numeric(colnames(mg))] = mean_glu
vals_nonsig = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = vals)

purple_palette = colorRampPalette(c("#FFFFFF", "purple"))(100)

gludist_nonsig = plot.hex.net.gg(cat = vals_nonsig$Value*10,
                                 altcols = purple_palette,
                                 are.unoccupied = which(vals_nonsig$Value == 0),
                                 remake = TRUE) +
  annotate("text",
           x = min(atr[which(vals_nonsig$Value != 0), 10]), # Centered hrzntl
           y = max(atr[which(vals_nonsig$Value != 0), 11]), # Near the top
           label = "I",
           size = 14/.pt,
           fontface = "bold",
           hjust = 0.5) +
  annotate("text",
           x = 0, # min(atr[which(vals_nonsig != 0), 10]),
           y = min(atr[which(vals_nonsig$Value != 0), 10])*0.55,
           label = paste0("No signaling: \n",
                          sum(vals_nonsig$Value > 0 &
                                vals_nonsig$Value < 0.001),
                          " starved cells"),
           color = "black", size = 2.7, hjust = 0.5)

gludist_sig = plot.hex.net.gg(cat = vals_stabsig$Value*10,
                altcols = purple_palette,
                are.unoccupied = which(vals_stabsig$Value == 0),
                remake = TRUE) +
  annotate("text",
           x = min(atr[which(vals_stabsig$Value != 0), 10]), # Centered hrzntl
           y = max(atr[which(vals_stabsig$Value != 0), 11]), # Near the top
           label = "H",
           size = 14/.pt,
           fontface = "bold",
           hjust = 0.5) +
  annotate("text",
           x = 0, # min(atr[which(vals_stabsig != 0), 10]),
           y = min(atr[which(vals_stabsig$Value != 0), 10])*0.55,
           label = paste0("Stable Oscillations: \n",
                          sum(vals_stabsig$Value > 0 &
                                vals_stabsig$Value < 0.001),
                          " starved cells"),
           color = "black", size = 2.7, hjust = 0.5)

################################################################################
### Combine Everything #########################################################
################################################################################

# Save the grid of plots to a TIFF file
tiff("Figures/Fig6.tiff",
     width = 7.5, height = 7.5, units = "in", res = 300)

# Create a 3-row layout with custom heights
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 3, ncol = 3,
  heights = unit(c(2, 1, 1), "null")  # Top 2 rows have height 1, bottom row has height 2
)))

# Helper function to define viewports
define_viewport <- function(row, col) {
  viewport(layout.pos.row = row, layout.pos.col = col)
}

# Place plots in the grid
print(Trac_minsig / Glu_minsig, vp = define_viewport(1, 1))
print(Trac_uncsig / Glu_uncsig, vp = define_viewport(1, 3))
print(Trac_stabsig / Glu_stabsig, vp = define_viewport(1, 2))

print(gludist_sig, vp = define_viewport(2, 3))
print(gludist_nonsig, vp = define_viewport(3, 3))

print(pp2, vp = define_viewport(2:3, 1:2))
dev.off() # Close the TIFF device






################################################################################
### Figure S8 ##################################################################
################################################################################

s8_nonsig = plot.hex.net.gg(cat = vals_nonsig$Value*10,
                                 altcols = purple_palette,
                                 are.unoccupied = which(vals_nonsig$Value == 0),
                                 remake = TRUE) +
  annotate("text",
           x = min(atr[which(vals_nonsig$Value != 0), 10]), # Centered hrzntl
           y = max(atr[which(vals_nonsig$Value != 0), 11]), # Near the top
           label = "A",
           size = 14/.pt,
           fontface = "bold",
           hjust = 0.5)

s8_minsig = plot.hex.net.gg(cat = vals_minsig$Value*10,
                            altcols = purple_palette,
                            are.unoccupied = which(vals_minsig$Value == 0),
                            remake = TRUE) +
  annotate("text",
           x = min(atr[which(vals_minsig$Value != 0), 10]), # Centered hrzntl
           y = max(atr[which(vals_minsig$Value != 0), 11]), # Near the top
           label = "B",
           size = 14/.pt,
           fontface = "bold",
           hjust = 0.5)

s8_stabsig = plot.hex.net.gg(cat = vals_stabsig$Value*10,
                              altcols = purple_palette,
                              are.unoccupied = which(vals_stabsig$Value == 0),
                              remake = TRUE) +
  annotate("text",
           x = min(atr[which(vals_stabsig$Value != 0), 10]), # Centered hrzntl
           y = max(atr[which(vals_stabsig$Value != 0), 11]), # Near the top
           label = "C",
           size = 14/.pt,
           fontface = "bold",
           hjust = 0.5)


s8_uncsig = plot.hex.net.gg(cat = vals_uncsig$Value*10 + 7,
                                 altcols = purple_palette,
                                 are.unoccupied = which(vals_uncsig$Value == 0),
                                 remake = TRUE) +
  annotate("text",
           x = min(atr[which(vals_uncsig$Value != 0), 10]), # Centered hrzntl
           y = max(atr[which(vals_uncsig$Value != 0), 11]), # Near the top
           label = "D",
           size = 14/.pt,
           fontface = "bold",
           hjust = 0.5)


combined = (s8_nonsig / s8_stabsig) | (s8_minsig / s8_uncsig)

ggsave("Figures/FigS7.tiff", combined,
       width = 5.2, height = 5.2, units = "in", dpi = 600)

### Calculate Gini coefficient

library(ineq)

Gini(vals_nonsig$Value)
Gini(vals_stabsig$Value)

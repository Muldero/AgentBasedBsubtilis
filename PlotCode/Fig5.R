# Get list of all RDS files in the folder
files = list.files("Revisions/Plot_Data_3500_30_slow/", pattern = "\\.RDS$",
                   full.names = TRUE)

# files = files[2:length(files)]

k_mat = readRDS("Revisions/K_e_Diffusion_Traj/K_e_diffusion_trajectory_3500.RDS")

starts = rev(seq((length(k_mat[1, ]) - 450), 151, by = -33))

(starts - 162) / 33

### Calculate base line

k_dat1 = readRDS("Revisions/Plot_Data_3500_30_slow/K_Synch_Data_dist_857.RDS")
k_dat2 = readRDS("Revisions/Plot_Data_3500_30_slow20/K_Synch_Data_dist_857.RDS")
k_data = cbind(k_dat1, k_dat2[, 3:12])

traj_mins = get.only.max.peaks(k_data$k_traj[401:800],
                               mean(k_data$k_traj[401:800]), top = FALSE) + 10

def_range = c()
def_iqr = c()
def_var = c()

for (j in 1:length(traj_mins)) {
  locs = 1:20
  for (k in 1:20) {
    mins = get.only.max.peaks(k_data[1:400, k + 2],
                              mean(k_data[1:400, k + 2]), top = FALSE)
    locs[k] = mins[which.min(abs(mins - traj_mins[j]))]
    if (traj_mins[j] != traj_mins[which.min(abs(traj_mins - locs[k]))]) {
      locs[k] = NA
    }
  }
  def_range = c(def_range, max(locs, na.rm = TRUE) - min(locs, na.rm = TRUE))
  def_iqr = c(def_iqr, IQR(locs, na.rm = TRUE))
  def_var = c(def_var, var(locs, na.rm = TRUE))
}

def_range = mean(def_range[2:(length(def_range) - 1)])
def_iqr = mean(def_iqr[2:(length(def_iqr) - 1)])
def_var = mean(def_var[2:(length(def_var) - 1)])

### Calculate for oscillations

traj_mins = get.only.max.peaks(k_data$k_traj[401:1200],
                               mean(k_data$k_traj[401:1200]), top = FALSE) + 10

range_mat = matrix(ncol = length(files) + 1, nrow = length(traj_mins))
var_mat = matrix(ncol = length(files) + 1, nrow = length(traj_mins))
iqrange_mat = matrix(ncol = length(files) + 1, nrow = length(traj_mins))

iqrange_mat[, 1] = range_mat[, 1] = var_mat[, 1] = traj_mins

id_vec = files

# Loop over files
for (i in 1:length(files)) {
  file = files[i]

  # Extract 4-digit number from filename
  id_match = regmatches(file, regexpr("\\d{3,4}(?=\\.RDS$)", file, perl = TRUE))
  id = as.integer(id_match)
  id_vec[i] = id

  if (id >= 850) {
    # Read the file
    k_dat1 = readRDS(file)
    k_dat2 = readRDS(paste0(
      "Revisions/Plot_Data_3500_30_slow20/K_Synch_Data_dist_", id, ".RDS"))

    k_data = cbind(k_dat1, k_dat2[, 3:12])

    traj_mins = get.only.max.peaks(k_data$k_traj[401:1200],
                                   mean(k_data$k_traj[401:1200]),
                                   top = FALSE) + 10

    for (j in 1:length(traj_mins)) {
      locs = 1:20
      for (k in 1:20) {
        mins = get.only.max.peaks(k_data[401:1200, k + 2],
                                  mean(k_data[401:1200, k + 2]), top = FALSE)
        locs[k] = mins[which.min(abs(mins - traj_mins[j]))]
        if (traj_mins[j] != traj_mins[which.min(abs(traj_mins - locs[k]))]) {
          locs[k] = NA
        }
      }
      range_mat[j, (id - 162)/33 + 2] = (max(locs, na.rm = TRUE) -
                                           min(locs, na.rm = TRUE))
      iqrange_mat[j, (id - 162)/33 + 2] = IQR(locs, na.rm = TRUE)
      var_mat[j, (id - 162)/33 + 2] = var(locs, na.rm = TRUE)
    }
  }
}

range_mat[, 2:length(range_mat[1, ])] =
  range_mat[, 2:length(range_mat[1, ])] / def_range
iqrange_mat[, 2:length(iqrange_mat[1, ])] =
  iqrange_mat[, 2:length(iqrange_mat[1, ])] / def_iqr
var_mat[, 2:length(range_mat[1, ])] =
  var_mat[, 2:length(range_mat[1, ])] #/ def_var

range_mat[which(is.infinite(range_mat), arr.ind = TRUE)] = NA
iqrange_mat[which(is.infinite(iqrange_mat), arr.ind = TRUE)] = NA

# range_mat = range_mat[, -(2:33)]

# ordering = order(as.numeric(id_vec))

# range_mat = cbind(range_mat[, 1, drop = FALSE], range_mat[, ordering + 1])
# var_mat = cbind(var_mat[, 1, drop = FALSE], var_mat[, ordering + 1])

plot(x = iqrange_mat[, 1], y = iqrange_mat[, 23])
# plot(x = var_mat[, 1], y = var_mat[, 2])


plot(y = colMeans(range_mat[12:15, 23:length(range_mat[1, ])]),
     x = sort(as.numeric(id_vec[as.numeric(id_vec) >= 850])), type = "l")
plot(y = colMeans(var_mat[12:15, 2:length(var_mat[1, ])]),
     x = sort(as.numeric(id_vec)))

k_ggdat = data.frame(
  Distance = sort(as.numeric(id_vec)),
  Effect = colMeans(iqrange_mat[12:15, 2:length(iqrange_mat[1, ])])
)

k_falloff_plot = k_ggdat[22:40, ] %>% ggplot(aes(x = Distance * 3 - 450)) +
  geom_point(aes(y = Effect), color = "gray", alpha = 0.5) +
  geom_smooth(aes(y = Effect), color = "black", size = 0.8, se = FALSE) +
  #geom_vline(xintercept = 201, color = "slategray", size = 2, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 1), n.breaks = 3) +
  scale_x_continuous(expand = expansion(mult = 0.02), limits = ) +
  labs(title = bquote(bold(B)),
       x = "Distance (um)",
       y = "Scaled Range") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11),
    title = element_text(size = 12)
  )

# ggsave(paste0(return_folder, "K_Falloff_Plot.tiff"), k_falloff_plot,
#        width = 5.2, height = 3, units = "in", dpi = 300)

k_dat1 = readRDS("Revisions/Plot_Data_3500_30_slow/K_Synch_Data_dist_857.RDS")
k_dat2 = readRDS("Revisions/Plot_Data_3500_30_slow20/K_Synch_Data_dist_857.RDS")

k_data = cbind(k_dat1, k_dat2[, 3:12])

p2_osc = ggplot(k_data[201:1200, ], aes(x = time + 200)) +
  geom_line(aes(y = rescale(k_traj, c(k_data[, 3:5])), lty = "Scaled K+"),
            color = "gold", size = 0.8) +
  geom_line(aes(y = k_osc.1, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.2, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.3, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.14, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.5, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.6, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.7, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.8, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.18, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = k_osc.10, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.5, size = 0.4) +
  scale_linetype_manual(values = c("Signaling Fraction" = "solid",
                                   "Scaled K+" = "21")) +
  scale_x_continuous(n.breaks = 6, expand = expansion(mult = 0.02)) +
  scale_y_continuous(n.breaks = 4) +
  labs(title = bquote(bold(A)),
       x = "Time",
       y = "Signaling Fraction",
       linetype = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11),
    title = element_text(size = 12),
    legend.position = c(0.9, .125),
    legend.margin = margin(0.1, 0, 0, 0.1, unit = "cm"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm")
  )

combo = p2_osc / k_falloff_plot + plot_layout(heights = c(3, 2))

combo

ggsave(paste0(return_folder, "Fig5.tiff"), combo,
       width = 7.5, height = 5, units = "in", dpi = 300)

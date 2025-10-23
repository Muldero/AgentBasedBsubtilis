# Get list of all RDS files in the folder
files = list.files("Revisions/Glu_Synch/",
                   pattern = "\\.RDS$", full.names = TRUE)

files = files[1:length(files)]

### Calculate base line

glu_data = readRDS(files[1])
traj_mins = get.only.max.peaks(glu_data$glu_traj[401:800],
                               mean(glu_data$glu_traj[401:800]), top = FALSE)

def_range = c()
def_iqr = c()
def_var = c()

for (j in 1:length(traj_mins)) {
  locs = 1:20
  for (k in 1:20) {
    mins = get.only.max.peaks(glu_data[1:400, k + 2],
                              mean(glu_data[1:400, k + 2]), top = FALSE)
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

traj_mins = get.only.max.peaks(glu_data$glu_traj[401:1200],
                               mean(glu_data$glu_traj[401:1200]), top = FALSE) - 15

range_mat = matrix(ncol = length(files) + 1, nrow = length(traj_mins))
iqr_mat = matrix(ncol = length(files) + 1, nrow = length(traj_mins))
var_mat = matrix(ncol = length(files) + 1, nrow = length(traj_mins))

iqr_mat[, 1] = range_mat[, 1] = var_mat[, 1] = traj_mins

id_vec = files

# Loop over files
for (i in 1:length(files)) {
  file = files[i]

  # Extract 4-digit number from filename
  id = as.numeric(sub(".*_(\\d+\\.?\\d*)\\.RDS$", "\\1", file))
  id_vec[i] = id

  # Read the file
  glu_data = readRDS(file)

  traj_mins = get.only.max.peaks(glu_data$glu_traj[401:1200],
                                 mean(glu_data$glu_traj[401:1200]), top = FALSE) - 15

  for (j in 1:length(traj_mins)) {
    locs = 1:20
    for (k in 1:20) {
      mins = get.only.max.peaks(glu_data[401:1200, k + 2],
                                mean(glu_data[401:1200, k + 2]), top = FALSE)
      locs[k] = mins[which.min(abs(mins - traj_mins[j]))]
      if (traj_mins[j] != traj_mins[which.min(abs(traj_mins - locs[k]))]) {
        locs[k] = NA
      }
    }
    range_mat[j, i + 1] = (max(locs, na.rm = TRUE) - min(locs, na.rm = TRUE))
    iqr_mat[j, i + 1] = IQR(locs, na.rm = TRUE)
    var_mat[j, i + 1] = var(locs, na.rm = TRUE)
  }
}

range_mat[, 2:length(range_mat[1, ])] =
  range_mat[, 2:length(range_mat[1, ])] / def_range
iqr_mat[, 2:length(range_mat[1, ])] =
  iqr_mat[, 2:length(range_mat[1, ])] / def_iqr
var_mat[, 2:length(range_mat[1, ])] =
  var_mat[, 2:length(range_mat[1, ])] / def_var

range_mat[which(is.infinite(range_mat), arr.ind = TRUE)] = NA

ordering = order(as.numeric(id_vec))

# iqr_mat = cbind(iqr_mat[, 1, drop = FALSE], iqr_mat[, ordering + 1])

#
# range_mat = cbind(range_mat[, 1, drop = FALSE], range_mat[, ordering + 1])
# var_mat = cbind(var_mat[, 1, drop = FALSE], var_mat[, ordering + 1])
#
# plot(x = range_mat[, 1], y = range_mat[, 20])
# # plot(x = var_mat[, 1], y = var_mat[, 2])
#
#
# plot(y = colMeans(range_mat[13:16, 2:length(range_mat[1, ])]),
#      x = sort(as.numeric(id_vec)), type = "l")
# plot(y = colMeans(var_mat[13:16, 2:length(var_mat[1, ])]),
#      x = sort(as.numeric(id_vec)))

glu_ggdat = data.frame(
  GluMag = as.numeric(id_vec), #sort(as.numeric(id_vec)),
  Effect = colMeans(iqr_mat[13:16, 2:length(iqr_mat[1, ])])
)

glu_falloff_plot = glu_ggdat[3:18, ] %>% ggplot(aes(x = GluMag)) +
  geom_point(aes(y = Effect), color = "black", alpha = 1) +
  # geom_smooth(aes(y = Effect), color = "black", size = 0.8, se = FALSE) +
  #geom_vline(xintercept = 201, color = "slategray", size = 2, alpha = 0.5) +
  scale_y_continuous(n.breaks = 3) +
  scale_x_log10(expand = expansion(mult = 0.02),
                breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10")) +
  labs(title = element_blank(), # bquote(bold(B)),
       x = "Glutamate Oscillation Magnitude (mM)",
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

ggsave(paste0(return_folder, "FigS6_quartile.tiff"), glu_falloff_plot,
       width = 5.2, height = 3, units = "in", dpi = 300)

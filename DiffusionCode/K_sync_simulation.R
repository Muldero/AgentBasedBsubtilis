age = 3000

k_mat = readRDS(
  "Plot_Code/Plot_Data/K_e_diffusion_trajectory_10000.RDS")

starts = rev(seq((length(k_mat[1, ]) - 450), 151, by = -33))

starts = starts[starts < 1600]

trajs = matrix(ncol = length(starts), nrow = 200)

k_traj_mat = matrix(ncol = length(starts), nrow = 800)

for (i in 1:length(starts)) {
  for (j in 801:1000) {
    trajs[j - 800, i] = mean(k_mat[j, starts[i]:(starts[i] + 300)])
  }
  trajs[, i] = trajs[, i] - min(trajs[, i])

  mins = get.only.max.peaks(trajs[, i], midline = mean(trajs[, i]), top = FALSE)
  k_traj = rep(trajs[mins[2]:(mins[length(mins)]), i], 20)

  k_traj_mat[, i] = k_traj[1:800] + 8
}

k_plot_data = data.frame(
  time = (-400):800,
  k_traj = c(rep(8, 401), k_traj_mat[, run_number])
)

for (repeater_var_10 in 1:10) {
  is.growing = FALSE
  atr = readRDS(paste0(read_folder, "POST_ATR/Post_atr_",
                       repeater_var_10, ".RDS"))
  atr = as.matrix(atr) # some of my code for ATRs doesn't work on dataframes
  trackers = readRDS(paste0(read_folder, "POST_TRACK/Post_track_",
                            repeater_var_10, ".RDS"))
  tht_mat = readRDS(paste0(read_folder, "THT/Tht_mat_",
                           repeater_var_10, ".RDS"))
  tht_vec = tht_mat[length(tht_mat[, 1]), ]
  tht_mat = rbind(tht_mat, matrix(0, nrow = 800, ncol = length(tht_mat[1, ])))

  age = length(trackers[, 1])
  unoccupied = which(atr[, 2] == -1)
  occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
  glu_radius = find.glu.radius()

  for (stepwise_timer in (age + 1):(age + 800)) {

    K_m = k_traj[stepwise_timer - age]

    atr[, 6:7] = glu.update(atr)[, 6:7]

    atr2 = atr
    atr2[, 4:5] = pot.update(atr)[, 4:5]
    atr2[, 8] = mem.potential.update(atr)[, 8]
    atr = atr2

    tht_vec = tht.update(tht_vec, atr)

    tht_vec[unoccupied] = -1

    tht_mat[stepwise_timer, ] = tht_vec

    atr = atr2

    # Correct values in unoccupied cells
    atr[unoccupied, c(5, 7, 8)] = 0
    atr[unoccupied, 2:3] = -1

    trackers = update.trackers(stepwise_timer)
  }

  print(repeater_var_10)

  post_traj_k_osc_1 = trackers$Out_tht_sig[(age - 400):(age + 800)]

  k_plot_data = cbind(k_plot_data, post_traj_k_osc_1)
  colnames(k_plot_data)[length(k_plot_data[1, ])] =
    paste0("k_osc.", repeater_var_10)
}

saveRDS(k_plot_data, paste0(plotdata_folder,
                            "K_Synch_Data_dist_", starts[run_number], ".RDS"))


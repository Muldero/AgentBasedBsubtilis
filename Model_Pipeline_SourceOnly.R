### Let's do this! #############################################################

################################################################################

# Initialize the biofilm

parent_matrix = cbind(1:7, c(NA, 1, 1, 1, 1, 1, 1)) # parents for cells in atr
colnames(parent_matrix) = c("Offspring", "Parent")

if (grow.all.the.way == TRUE) {
  atr = grow.a.biofilm2(method = "entire")
} else {
  atr = grow.a.biofilm2(method = "first_layer") # make.hex.atr()
  for (i in 1:ceiling(init_radius/gens_per_tick)) {
    atr = grow.a.biofilm2(method = "add_layer")
  }
}

atr_saver = atr
atr = atr_saver

################################################################################

unoccupied = which(atr[, 2] == -1)
occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
glu_radius = find.glu.radius()

# record potassium, glutamate and memebrane potential for 100 random cells
pot_traj = matrix(nrow = 100, ncol = 0)
glu_traj = matrix(nrow = 100, ncol = 0)
mem_pot_traj = matrix(nrow = 100, ncol = 0)

traj.created = FALSE
is.growing = TRUE

tht_mat = matrix(0, nrow = ticks, ncol = length(atr[, 1]))
tht_vec = rep(0, length(atr[, 1]))

sigs_mat = matrix(0, nrow = ticks, ncol = length(atr[, 1]))

# Now for the actual simulation
stepwise_timer = 0
for (stepwise_timer in (stepwise_timer + 1):ticks) {

  # Do we do a potassium shock?
  if (stepwise_timer >= min(K_shock) & stepwise_timer < max(K_shock + 5)) {
    K_m = kshock_level
  } else {
    K_m = 8
  }

  # Simultaneously update potassium, glutamate, and membrane potential
  atr[, 6:7] = glu.update(atr)[, 6:7]

  for (subtick in 1) {
    atr2 = atr
    atr2[, 4:5] = pot.update(atr, tickpper = tpp * 1)[, 4:5]
    # atr2[, 6:7] = glu.update(atr)[, 6:7]
    atr2[, 8] = mem.potential.update(atr, tickpper = tpp * 11)[, 8]
    atr = atr2
  }

  tht_vec = tht.update(tht_vec, atr)

  tht_vec[unoccupied] = -1

  tht_mat[stepwise_timer, ] = tht_vec

  # Correct values in unoccupied cells
  atr[unoccupied, c(5, 7, 8)] = 0
  atr[unoccupied, 2:3] = -1


  # While growing
  if (length(occupied)/length(atr[, 1]) < prop_occupied) {
    # growth the biofilm by 1/40 of a layer
    atr = grow.a.biofilm2(method = "add_layer")

    # update occupied cells and inner/outer boundary
    unoccupied = which(atr[, 2] == -1)
    occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
    glu_radius = find.glu.radius()

  } else if (is.growing == TRUE) {
    # first tick after growth stops perform some stuff

    is.growing = FALSE # stop growth

    cells_to_track = c(sample(find.external(), 800),
                       sample(find.internal(), 200))
    # all of the above are just various recording mechanisms

    if (track.min.glu == TRUE) {
      # make minimum glutamate tracking matrix
      min_glu = atr[occupied, 7]
    }
  } else {
    # After growth starts

    if (track.min.glu == TRUE) {
      # update minimum glutamate tracker
      min_glu = rbind(min_glu, atr[occupied, 7])
    }

    # Update single cell traces
    pot_traj = cbind(pot_traj, atr[cells_to_track[1:100], 5])
    glu_traj = cbind(glu_traj, atr[cells_to_track[1:100], 7])
    mem_pot_traj = cbind(mem_pot_traj, atr[cells_to_track[1:100], 8])
  }

  # update trackers (more diagnostics and output)
  trackers = update.trackers(stepwise_timer)
  plot.singlecell.trajs(n = 32, do.print = FALSE)

  print(stepwise_timer)

  if (stepwise_timer %% 7 == 0) {
    par(mfrow = c(4, 1))
    # plot(density(tht_vec[occupied]), xlab = "Scaled ThT", main = stepwise_timer)

    plot(trackers$Tot_K_e[1:stepwise_timer], type = "l", xlab = "Time",
         ylab = "K_e", main = "")

    plot(trackers$Out_tht_sig[1:stepwise_timer], type = "l", xlab = "Time",
         ylab = "ThT Signalers", main = "")
    abline(h = 0.43)

    plot(density(Eff(atr[occupied, 8])), xlab = "F(V)", main = "")

    plot(trackers$Out_fv_sig[1:stepwise_timer], type = "l", xlab = "Time",
         ylab = "F_V Signalers", main = "")
    abline(h = 0.43)

    # plot(density(atr[occupied, 5]), xlab = "K_i", main = "")
  }
}

if (track.min.glu == TRUE & exists("min_glu")) {
  min_glu = rbind(occupied, min_glu)
}






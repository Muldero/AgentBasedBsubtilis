flux = function(D = "diffusion coefficient",
                c1 = "concentration 1",
                c2 = "concentration 2",
                delta = "distance over which diffusion takes place",
                a = "Area (of the border)") {
  J = D * (c1 - c2) / delta # mM/um^2/hr
  J = 10^(-3) * J # pico mols / um^2 / hr
  J = J * time_step # pico moles / um^2 / time_step
  return(J*a) # total pico moles in a tick
}

# dist_vec = seq(2000, 10000, by = 500)
#
# distance = dist_vec[run_number]

distance = 10000

edge = 150 + ceiling(distance/3)

total = edge + 450 # give an extra 150 before setting to 8

################################################################################
### Realist K_e version ########################################################
################################################################################

trackers = readRDS("Output/POST_TRACK/Post_track_1.RDS")

k_e = trackers$Tot_K_e[2001:3000]

k_e_adjust = k_e # ((k_e - min(k_e)) / (max(k_e) - min(k_e))) * 7 + 5

plot(k_e_adjust, type = "l")

# each tick is 1.2 minutes, total time for 300 oscillations is 6 hours
# hypothetical second biofilm is at cell 817 (border)

### Set up matrix

K_m = 8

D_e = (7/4)*10^6 # um^2/hr
D_o = 7*10^6

diff_mat = data.frame(
  "Depth" = 0:(total - 1),
  "Volume" = (1:total) * 3, # assuming cells 1x1x3 um, units um^3
  "Concentration" = c(rep(k_e_adjust[1], 150), rep(K_m, total - 150)), # units mMol
  "picomoles" = c(rep(k_e_adjust[1], 150), rep(K_m, total - 150)) *
    (1:total) * 3 * 10^(-3), # picomoles
  "D" = c(rep(D_e, 150), rep(D_o, total - 150))
)

time_step = 10^-8 # ticks per hour

### Set up records
record_per_step_mat = matrix(nrow = 2*10^6, ncol = total) # time-steps per tick
tick_record_mat = matrix(nrow = 1000, ncol = total)

### Run

for (i in 1:(1000 * 2 * 10^6)) {
  tick_count = floor(i/(2 * 10^6)) + 1

  # diff_mat_temp = diff_mat

  diff_vec = flux(
    (diff_mat$D[1:(total - 1)] + diff_mat$D[2:total]) / 2, # diffusion coefficient
    diff_mat$Concentration[1:(total - 1)], # concentration 1
    diff_mat$Concentration[2:total], # concentration 2
    rep(3, total - 1), # distance, each layer 3 um thick
    1 # diff_mat$Volume[1:(total - 1)] # um, should this be in units of um^2 or um^3??
  ) # flux in picomoles per one unit border

  diff_vec = diff_vec * (1:(total - 1) + 0.5)

  diff_mat$picomoles[1:(total - 1)] = diff_mat$picomoles[1:(total - 1)] -
    diff_vec
  diff_mat$picomoles[2:total] = diff_mat$picomoles[2:total] + diff_vec

  diff_mat$Concentration = 10^3 * diff_mat$picomoles / diff_mat$Volume

  diff_mat$Concentration[total] = K_m
  diff_mat$picomoles[total] = K_m * 10^(-3) * total * 3

  diff_mat$Concentration[1:150] = k_e_adjust[tick_count]
  diff_mat$picomoles[1:150] = k_e_adjust[tick_count] * 10^(-3) * 1:150 * 3

  record_per_step_mat[((i - 1) %% (2 * 10^6)) + 1, ] = diff_mat$Concentration

  if (i %% (2 * 10^6) == 0) {
    tick_record_mat[tick_count - 1, ] = colMeans(record_per_step_mat)
    print(paste0("Tick:", tick_count - 1, " | Time: ", Sys.time()))
  }
}

### Save Output
saveRDS(tick_record_mat,
        "PlotCode/Plot_Data/K_e_diffusion_trajectory_10000.RDS")






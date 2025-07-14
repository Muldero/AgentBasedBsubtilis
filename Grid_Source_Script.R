### Library ####################################################################

# setwd("~/Desktop/Post_Reviewer/EndeavourScripts")

# Error handling
options(error = function() {
  cat("\n❌ ERROR DETECTED — dumping diagnostic data to ERRORS/\n")
  dump.frames("ERRORS/error_dump", to.file = TRUE)
  save.image("ERRORS/error_workspace.RData")
  q("no", status = 1)
})

source("Packages.R", verbose = FALSE)

sig_cutoff = 0.035
internal_adjust = 8
V_0T = -95

options(max.print = 3000, digits = 8)
options(scipen = 8)

plot.gif = FALSE # plot a gif of the model
grow.all.the.way = FALSE # grow the entire biofilm (remove the growth phase)

# memory intensive (!) diagnostic trackers
track.min.glu = TRUE
track.all.sigs = FALSE

ticks = 3000 # how many ticks to run the model for

tpp = 100 # ticks per period (of 2 hours) so divide by tpp/2 for time-dependence

gens_per_tick = 1/40 # number of generations (amount of growth) per tick

source("Functions.R", verbose = FALSE)

return_folder = paste0("Grid_Output/")

### Variables ##################################################################

r = 150 # radius of network
init_radius = 40
prop_occupied = 0.75 # proportion of the network that should be occupied if grow

# Observed (true) parameters
K_shock = K_shock_default = -10
K_m = 8 # basal potassium concentration in media
G_m = 30 # basal glutamate concentration in media
V_0 = -86 # base membrane potential

# Potassium parameters
F_param = 0.05/0.23 # default 0.05
D_p = 0.12*2.35 # 1.5 # default 0.12

G_m_vec = rep(G_m, r*6 + 6)

# Parameters (drawn from Martinez-Corral (2019))
alpha_g = 24 # glu uptake constant
delta_g = 5.7 # 4.8 # glu degradation constant
g_k = 60 # 70 # default 70
g_l = 18 # 18 in MC< 12 in Prindle

# Electrical parameters
V_K0 = 25.8 # Nernst potential prefactor
V_L0 = -93.5 # Basal leak potential - https://www.genomicscience.energy.gov/abstra
# ct/fluorescence-lifetime-based-imaging-o
# f-bacillus-subtilis-membrane-potential/
V_0 = -86 # base membrane potential: V_L0 - V_0 >> 0 harder to be effected by
# depolarization. V_L0 - V_0 = 0 easy to
# be effected by depolarization
d_l = 7.6 # default 4 (MC) or 8 (P + F)

k_method = "MC" # "MC" "PF"

# Threshold Parameters
upper_thresh = upper_thresh_default = 3
lower_thresh = lower_thresh_default = 0
thresh_herit = 1
sd2 = sd2_default = 1 # standard deviation of genotype

# source("Master_Model_Pipeline_SourceOnly.R")

args = commandArgs(trailingOnly = TRUE)
run_count = as.numeric(args[1])

# Define the vector
min_values = seq(-0.6, 0.6, by = 0.1)
max_values = seq(2.5, 3.5, by = 0.1)

# Create combinations of upper and lower where upper >= lower
thresh_mat = expand.grid(lower = min_values, upper = max_values)
thresh_mat = thresh_mat[thresh_mat$upper >= thresh_mat$lower, ]

# Add the sigma column
thresh_mat$sigma = ifelse(thresh_mat$upper == thresh_mat$lower, 0,
                          (thresh_mat$upper - thresh_mat$lower) / 3)


upper_thresh = thresh_mat$upper[run_count]
lower_thresh = thresh_mat$lower[run_count]
sd2 = thresh_mat$sigma[run_count]

source("Master_Model_Pipeline_SourceOnly.R")

# Record Output
saveRDS(atr, paste0(return_folder, "POST_ATR/Post_atr_", run_count, ".RDS"))

saveRDS(trackers, paste0(return_folder, "POST_TRACK/Post_track_",
                         run_count, ".RDS"))

# saveRDS(tht_mat, paste0(return_folder, "THT/Tht_mat_",
#                         run_count, ".RDS"))
#
# saveRDS(sigs_mat, paste0(return_folder, "SIG_MAT/Sigs_mat_",
#                          run_count, ".RDS"))

saveRDS(pot_traj, paste0(return_folder, "Trajectories/Potassium_Trajectories",
                         run_count, ".RDS"))

saveRDS(glu_traj, paste0(return_folder,
                         "Trajectories/Glutamate_Trajectories_",
                         run_count, ".RDS"))

saveRDS(mem_pot_traj,
        paste0(return_folder, "Trajectories/Membrane_Potential_Trajectories_",
               run_count, ".RDS"))

saveRDS(find.external(), paste0(return_folder, "Other/External_cells_",
                                run_count, ".RDS"))

saveRDS(min_glu, paste0(return_folder, "MIN_GLU/Minimum_Glutamate_",
                        run_count, ".RDS"))

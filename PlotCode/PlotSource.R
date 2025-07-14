### Library ####################################################################

setwd("~/Desktop/Correct_V0")

source("Packages.R", verbose = FALSE)

sig_cutoff = 0.035
internal_adjust = 8

options(max.print = 3000, digits = 8)
options(scipen = 8)
plot.gif = FALSE
grow.all.the.way = FALSE
track.min.glu = FALSE
track.all.sigs = FALSE
ticks = 3000
tpp = 100
gens_per_tick = 1/40

source("Functions.R", verbose = FALSE)
read_folder = paste0("Output/")
plotdata_folder = paste0("PlotCode/Plot_Data/")
return_folder = paste0("Figures/")

parent_matrix = cbind(1:7, c(NA, 1, 1, 1, 1, 1, 1))

### Variables ##################################################################

r = 150
init_radius = 40
prop_occupied = 0.75
K_shock = K_shock_default = -10
K_m = 8
G_m = 30
F_param = 0.05/0.23
D_p = 0.12*2.35
G_m_vec = rep(G_m, r*6 + 6)
alpha_g = 24
delta_g = 5.7
g_k = 60
g_l = 18
V_K0 = 25.8
V_L0 = -93.5
V_0 = -86
V_0T = -95
d_l = 7.6
k_method = "MC"
upper_thresh = upper_thresh_default = 3
lower_thresh = lower_thresh_default = 0
thresh_herit = 1
sd2 = sd2_default = 1 # standard deviation of genotype

args = commandArgs(trailingOnly = TRUE)
run_number =  as.numeric(args[1])

################################################################################

source("PlotCode/RelevantPlot.R")

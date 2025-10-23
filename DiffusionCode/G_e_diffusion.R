# We are eventually going to have 4 models, of increasing complexity. First,
# we'll simulate a slice of the biofilm (area of each layer is constant), with
# constant glutamate drain and no potassium or anything like that.
# Second, we'll do the same, but simulate an arc instead (all that changes is
# the area of each layer increases as we go out).
# Next, if necessary, we'll reproduce both of those, but simulating glutamate,
# potassium, and membrane potential as in Martinez-Corral 2019

setwd("~/Desktop/Correct_V0/Revisions")

### Functions ##################################################################

# Function to calculate flow of a molecule from one layer of the simulation to
# another. Takes the diffusion coefficient, concentrations in each layer, the
# depth of each layer, and the area of the border. Returns picomoles, which
# must be converted appropriately.
flux = function(D = "diffusion coefficient", # um^2/hr
                c1 = "concentration 1", # mM
                c2 = "concentration 2", # mM
                delta = "distance over which diffusion takes place", # um
                l = "Length (of the border)") { # um
  J = D * (c1 - c2) / delta # mM/um/hr
  J = l * 10^(-6) * J # * delta # pmol / hr         ######################################## mM to picomol needs volume not just 10-6 right?
  J = J * time_step # pmol / time_step
  return(J) # total picomoles crossing in a tick
}

### Diffusion parameters #######################################################

# K+ coefficient citation: https://www.aqion.de/site/diffusion-coefficients
# Glu coefficient citation: https://pmc.ncbi.nlm.nih.gov/articles/PMC3133640/
# Biofilm coefficient scaling: Stewart 2003 Diffusion in Biofilms

# Dk_m = 7*10^6 # um^2/hr - potassium diffusion coefficient in the media
# Dk_b = Dk_m/4 # um^2/hr - potassium diffusion coefficient in the biofilm

Dg_m = 2.7*10^6 # um^2/hr - glutamate diffusion coefficient in the media
Dg_b = Dg_m/4 # um^2/hr - glutamate diffusion coefficient in the biofilm

G_m = 30 # mM

### Other parameters ###########################################################

layer_depth = 3

width = 131 * (3 / layer_depth)

total = width + 1 * (3 / layer_depth)

### Time

# in original model, each tick was 1.2 minutes
# for accurate diffusion, we need ticks to be 10^-8/hr
# replace tpp/2 w time_step in equations

time_step = 10^-8 # ticks/hr

# Parameters (drawn from Martinez-Corral (2019))
alpha_g = 20 # 24 # glu uptake constant
delta_g = 2.69 # 3.25 #4.8 # 4.8 # glu degradation constant

g_degrade_rate = delta_g*time_step

### Matrix setup ###############################################################

# Now we add in scaling
# Inner 'width' layers are biofilm (numbering from 1:width), then a single outer
# layer that we reset the values for every tick.

diff_mat = data.frame(
  "Layer" = 1:total,
  "Volume" = (1:total) * layer_depth, # assuming cells 1x1xlayer_depth um, units um^3
  "G_In_Concentration" = rep(0, total), # units mMol
  "G_In_pmols" = rep(0, total), # picomoles
  "G_Ex_Concentration" = c(rep(1, width), rep(G_m, total - width)), # units mM
  "G_Ex_pmols" = c(rep(1, width), rep(G_m, total - width)) *
    layer_depth * 10^(-6) * 1:total, # picomoles
  "Dg" = c(rep(Dg_b, width), rep(Dg_m, total - width))
)

### Simulations ################################################################

### Set up records
record_per_step_ex_mat = matrix(nrow = 2*10^5, ncol = total) # time-steps per tick
tick_record_ex_mat = matrix(nrow = 10000, ncol = total)

record_per_step_in_mat = matrix(nrow = 2*10^5, ncol = total) # time-steps per tick
tick_record_in_mat = matrix(nrow = 10000, ncol = total)

### Run

par(mfrow = c(2, 1))
for (i in 1:(10 * 2 * 10^6)) {
  tick_count = floor(i/(2 * 10^5)) + 1

  # Handle diffusion
  diff_vec = flux(
    (diff_mat$D[2:total] + diff_mat$D[1:(total - 1)]) / 2, # diff coefficient
    diff_mat$G_Ex_Concentration[2:total], # concentration 1
    diff_mat$G_Ex_Concentration[1:(total - 1)], # concentration 2
    rep(layer_depth, total - 1), # distance, each layer 3 um thick
    ((diff_mat$Volume[2:total] + diff_mat$Volume[1:(total - 1)]) / 2) / layer_depth # area, not volume here
  ) # flux in picomoles per one unit border

  diff_mat$G_Ex_pmols[2:total] = diff_mat$G_Ex_pmols[2:total] - diff_vec
  diff_mat$G_Ex_pmols[1:(total - 1)] =
    diff_mat$G_Ex_pmols[1:(total - 1)] + diff_vec
  diff_mat$G_Ex_pmols[total] = G_m * 10^(-6) * diff_mat$Volume[total]
  diff_mat$G_Ex_Concentration = 10^6 * diff_mat$G_Ex_pmols / diff_mat$Volume

  # Now do uptake and absorption simultaneously (ignore vol difference for now)

  absorption_vec = alpha_g * time_step * diff_mat$G_Ex_Concentration[1:width] /
    (0.75 + diff_mat$G_Ex_Concentration[1:width])
  degredation_vec = diff_mat$G_In_pmols[1:width] * g_degrade_rate * 100

  av_pmols = absorption_vec * diff_mat$Volume[1:width] * 10^(-6) * 100 # correction factor

  diff_mat$G_Ex_pmols[1:width] = diff_mat$G_Ex_pmols[1:width] - av_pmols
  diff_mat$G_In_pmols[1:width] = diff_mat$G_In_pmols[1:width] +
    av_pmols - degredation_vec

  diff_mat$G_In_pmols[diff_mat$G_In_pmols < 0] = 0
  diff_mat$G_Ex_pmols[diff_mat$G_Ex_pmols < 0] = 0

  diff_mat$G_In_Concentration = 10^6 * diff_mat$G_In_pmols / diff_mat$Volume
  diff_mat$G_Ex_Concentration = 10^6 * diff_mat$G_Ex_pmols / diff_mat$Volume

  record_per_step_ex_mat[((i - 1) %% (2 * 10^5)) + 1, ] =
    diff_mat$G_Ex_Concentration
  record_per_step_in_mat[((i - 1) %% (2 * 10^5)) + 1, ] =
    diff_mat$G_In_Concentration

  if (i %% (2 * 10^5) == 0) {
    tick_record_ex_mat[tick_count - 1, ] = colMeans(record_per_step_ex_mat, na.rm = TRUE)
    tick_record_in_mat[tick_count - 1, ] = colMeans(record_per_step_in_mat, na.rm = TRUE)

    plot(tick_record_ex_mat[tick_count - 1, ][1:(width)],
          xlab = "Distance from Center", ylab = "External |G|",
          main = i/(2 * 10^5))
    # abline(v = width)
    plot(y = tick_record_in_mat[tick_count - 1, ][1:(width)],
         x = (1:width),
         xlab = "Distance from Center", ylab = "Internal |G|",
         main = i/(2 * 10^5))
    # abline(v = width)
    print(paste0("Tick: ", tick_count, " | Time: ", Sys.time()))
  }
}
par(mfrow = c(1, 1))

### Save Output
saveRDS(tick_record_ex_mat,
        "ModelOutput/G_diffusion_external_trajectory.RDS")
saveRDS(tick_record_in_mat,
        "ModelOutput/G_diffusion_internal_trajectory.RDS")




################################

glm = readRDS("ModelOutput/GlutamateLayerMatrix.RDS")

tick_record_ex_mat = readRDS("ModelOutput/G_diffusion_external_trajectory.RDS")
tick_record_in_mat = readRDS("ModelOutput/G_diffusion_internal_trajectory.RDS")

tick_count = 100
width = 131

plot(glm[800, 1:131], main = "Model", ylab = "Internal |G|",
     xlab = "Distance from Center")

plot(tick_record_in_mat[tick_count - 1, ][1:(width)],
     xlab = "Distance from Center", ylab = "Internal |G|",
     main = "Diffusion Simulation")

library(ggplot2)
library(tidyverse)
library(patchwork)

p1 = data.frame(glm = glm[1000, 1:131],
                diff_dat = tick_record_in_mat[tick_count - 1, ][1:(width)],
                layer = 1:131) %>%
  ggplot() +
  geom_line(aes(y = glm, x = layer, lty = "Heuristic"),
            color = "purple") +
  geom_line(aes(y = diff_dat, x = layer, lty = "Explicit"),
            color = "purple") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.position = c(.87, .21),
    legend.spacing = unit(0, "cm"),
    legend.margin = margin(0, 0, 0, 00, unit = "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  labs(x = "Cell Layer",
       y = "Internal Glutamate (mM)",
       lty = "Diffusion \nModel")
p1

ggsave("SimComparison.tiff", p1, width = 5.2, height = 3,
       units = "in", dpi = 300)

scale_factor = max(tick_record_in_mat[tick_count - 1, ][1:(width)]) /
  max(tick_record_ex_mat[tick_count - 1, ][1:(width)])

p2 = data.frame(ex_glu = tick_record_ex_mat[tick_count - 1, ][1:(width)],
                in_glu = tick_record_in_mat[tick_count - 1, ][1:(width)],
                layer = 1:131) %>%
  ggplot() +
  geom_line(aes(y = ex_glu*scale_factor, x = layer),
            color = "#CC4F00", lty = "dashed", size = 0.5) +
  geom_line(aes(y = in_glu, x = layer),
            color = "purple", lty = "solid", size = 0.5) +
  scale_y_continuous(
    name = "Cellular Glutamate (mM)",
    sec.axis = sec_axis(trans = ~./scale_factor,
                        name = "Extracellular Glutamate (mM)")
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title.y.left = element_text(color = "purple"),
    axis.text.y.left = element_text(color = "purple"),
    axis.title.y.right = element_text(color = "#CC4F00"),
    axis.text.y.right = element_text(color = "#CC4F00"),
  ) +
  labs(x = "Cell Layer")
       # y = "Internal Glutamate (mM)",
       # title = bquote(bold("A")))

ggsave("SimGlu.tiff", p2, width = 5.2, height = 3,
       units = "in", dpi = 300)


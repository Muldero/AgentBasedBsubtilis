This repository contains the code used in the paper "An Agent-Based Model of
  Metabolic Signaling Oscillations in Bacillus subtilis Biofilms"

Authors: Obadiah Mulder, Maya Peters Kostman, Abdulrahmen Almodaimegh, Michael
Edge, and Joseph Larkin

For questions please contact Obadiah Mulder at omulder@usc.edu

Code was run in RStudio 2023.03.0+386 "Cherry Blossom" Release with
R version 4.3.2 (2023-10-31) -- "Eye Holes"

### License ####################################################################

This code is licensed under the [Creative Commons Attribution 4.0 International
  License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You are free to:
- **Share**: Copy and redistribute the material in any medium or format.
- **Adapt**: Remix, transform, and build upon the material for any purpose, even
  commercially.

Under the following terms:
- **Attribution**: You must give appropriate credit, provide a link to the
  license, and indicate if changes were made. You may do so in any reasonable
  manner, but not in any way that suggests the licensor endorses you or your use.

For more information, visit https://creativecommons.org/licenses/by/4.0/

### Before running #############################################################

This code requires a folder tree to store return data. There must be a directory
  "Output" in the main directory, and it has to contain "Other", "POST_ATR",
  "POST_TRACK", "SIG_MAT", and "THT". A second directory "Grid_Output" also 
  needs sub-directories "POST_ATR", "POST_TRACK", "Trajectories", "Other", and
  "MIN_GLU". The "PlotCode" directory needs a "Plot_Data" directory, with 
  "Glu_Sync" and K_"Sync" sub-directories. The main directory also needs a 
  "Figures" directory. I apologize for the inconvinience, GitHub does not 
  support empty folders.

RequiredPackages.csv contains a list of all the requisite pacakges with version
  numbers. Newest versions of packages are automatically installed and loaded,
  but if this causes errors then the appropriate version may be required. If
  not accessible through CRAN they are available on request.

### Files ######################################################################

SimulationScript.R - Running this file will load in all packages and run code to
  produce the data in our paper

Functions.R - Sourced in CompleteSourceScript.R. All of the functions required
  to run and analyze these simulations.

Model_SourceOnly.R - Sourced in CompleteSourceScript.R. This script contains the
  code to perform simulations. Must be run by sourcing
  SimulationScript.R.

Packages.R - Sourced in CompleteSourceScript.R. Loads all of the required
  packages. Sets library to the Packages directory.

Grid_Source_Script.R - This will run many simulations with different bounds on
  stress thresholds. These results are used in Figs 5, 6, and S7.

K_diffusion.R - Calculates the oscillations of potassium caused by a signaling
  biofilm over several minutes out to a distance of 10 mm.

K_sync_simulation.R - Given a distance, takes the return from K_diffusion.R and
  runs 10 simulations of a biofilm beign exposed to the K+ oscillations caused
  by a second biofilm signaling at the given distance. These results are used to
  create Figure 5. PlotSource.R must be run before this.

### Directories ################################################################

Output - This contains several empty directies. Where all of the data generated
  by running CompleteSourceScript.R is stored.

Packages - Contains all of the source code for the packages used in this code.

PlotCode - Code to generate the plots and other summary data reported in this
  paper. Scripts in that folder do not run independently. PlotSource.R must be
  run before any of the plotting code will work.

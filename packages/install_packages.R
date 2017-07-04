### Run this on first use of DMC to install requried packages. Run from RStudio
### in order to get help with dependencies.

# Standard packages from CRAN
install.packages("loo") # For WAIC and looaic calculation
install.packages("hypergeo") # For population plausible values
install.packages("statmod") # Wald model
install.packages("rtdists") # For standard model distirbution functions
install.packages("pracma")  # For gng and stop signal robust integration
install.packages("snowfall") # Parallel processing
install.packages("rlecuyer") # Parallel processing
install.packages("numDeriv") # Prior transformations
install.packages("vioplot") # Stop signal graphs
install.packages("ggplot2") # For fancy graphs
install.packages("gridExtra") # For fancy graphs

# Note: Before running commands scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 

# Modifed packages distributed with DMC.

## Note that on some systems you will need to have admin rights to install 
## the followingl packages sucessfully and/or have your PATH set correctly. 
## It is also sometimes necessary to launch RStudio from an admin account.

# This modificaiton of coda allows for plotting of priors with plot.dmc
install.packages("packages/coda_0.18-1.3.tar.gz",repos=NULL,type="source")
## -- Installation Note
## On RStudio, you may use the "Install" button in the "Packages" tab.
## From its drop-down menu, choose "Package Archive File (.tar.gz)" to
## browser the location you store the source file "coda_0.18-1.3.tar.gz"


## DEVELOPMENT STUFF NOT FOR GENERAL USE

# # If you are using PDA then you will need this package
# install.packages("packages/cpda_0.0.2.6.tar.gz",repos=NULL,type="source")
# 
# # On the command line of a linux system (the first command enables GCC recent 
# # version, may vary or not be needed depending on the system)
# # scl enable devtoolset-4 bash
# # R CMD INSTALL -l ~/R/x86_64-redhat-linux-gnu-library/3.2/ packages/cpda_0.0.2.6.tar.gz
# 
# # If you are using CircularDDM then you will need this package
# install.packages("dmc/packages/CircularDDM_0.1.5.tar.gz",repos=NULL,type="source")

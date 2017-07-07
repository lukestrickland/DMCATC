
# Run Model
rm(list=ls()) 
setwd("~/Modelling")
source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")
setwd("~/Modelling/x1/samples")

run.grid.dmc("E1.1sdv.samples",model.dir="LBA",model.file="lbaN_B.R",user="rjb779",n.add=40,wall.hours=30,GB=5)  
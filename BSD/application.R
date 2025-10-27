#-------------------------------------------------------------------
#-------------------------
# Application to real data: 

# Data from the Regional Labor Court of the 5th Region (Bahia, Brazil).
## This data frame contains the following columns:
### VARA: Regional Labor Court
### TC: Congestion Rate
### IC: Conciliation Index

#---------------------
# Load the functions
source("functions.R")

#--------------
library(gamlss)
library(readxl)
setwd('...') # Set the working directory.
dados <- read_excel('.../trt5.xlsx') #  Load the data.
head(dados,5)

y1 <- dados$TC
y2 <- dados$IC
y <- data.frame(y1,y2)

est <- estimBSB(y, dist = "Simplex", copula = "FGM")
summaryBSBC(est)

#--------------------------------------
# Calculation of the joint expectation:

theta <- c(est$parameters[1],est$parameters[2],est$parameters[3],
           est$parameters[4], est$parameters[5])

eBsimplex(theta, method = "struve")
eBsimplex(theta, method = "bessel")  
#----------------------------------
#
# Plot BSD Joint Density

param <- c(est$parameters[1],est$parameters[2],est$parameters[3],
           est$parameters[4], est$parameters[5])

plotBSD(y1, y2, param, dist = "Simplex",
          copula_type = "FGM",
          plot_type = "surface3d",
          title = NULL, n_grid = 150)

plotBSD(y1, y2, param, dist = "Simplex",
          copula_type = "FGM",
          plot_type = "contour",
          title = NULL, n_grid = 150)
#------------------------------------
#------------------------------------------------------------------

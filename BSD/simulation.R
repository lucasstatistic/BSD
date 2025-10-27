#--------------------
# Load the functions
#
source("functions.R")
#------------------
# Simulation Study:

sim <- rBSBC(150, copula = "FGM", theta = 0.5,  dist = "Simplex",
             par1 = c(0.5, 7.1), par2 = c(0.5, 5))

y <- data.frame(sim[,1],sim[,2])
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
#
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
#----------------------------------
#------------------------------------------------------------------
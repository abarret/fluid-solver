Solves the unsplit Navier-Stokes equations using the projection method as a preconditioner as demonstrated in:

Griffith, B. E. (2009). An accurate and efficient method for the incompressible Navier-Stokes equations using the projection method as a preconditioner. Journal of Computational Physics, 228(20), 7565–7595. https://doi.org/10.1016/j.jcp.2009.07.001


All discretizations are second order finite differences except the advection code. Currently, the advection u*grad(u) is discretized using a first order wave propagation method with WENO reconstructions.

Currently only uses periodic boundary conditions. The infrastructure is there for physical boundary conditions, but not yet implemented.

Main driver is found in main_driver.m

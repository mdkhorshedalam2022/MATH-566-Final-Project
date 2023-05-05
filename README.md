# MATH 566: Final Project

This program solves the 2D PFAS transport equation using MATLAB's built-in Gaussian Elimination Method and Generalized Minimal Residuals Method. The goal of this program is to simulate the transport of PFAS to understand the impact of different transport processes, including diffusion, advection, solid-phase adsorption, and air-water adsorption. The model simulates flow for a domain (length 2.1 m and thickness 2 m) of Accusand soil. The boundary walls of the domain are impermeable, but there is a midpoint inlet on the left boundary wall and three outlets on the right boundary wall.

The governing equation, boundary conditions, mathematical formulation, and implementation of the numerical scheme are discussed in the file attached named [Project Report.pdf](./Project Report.pdf).

## Configure Instructions
To run the program, one needs to download the files in the code folder, then open the [main.m](./main.m) in MATLAB interface. The program requires the user to provide values for certain parameters to suit the details of the desired analysis, though for a few parameters a default value can be chosen. 

## Output
The program will generate two figures. The first figure will show the contour lines of PFAS concentration; the second figure will show the plot of the time history of PFAS concentration at three specific nodes. Also, this program will print the elapsed time based on the chosen numerical method. 

In addition, to generate the computational time comparison plot of the Gaussian elimination method and Generalized minimal residual method, one needs to download and run [comparisonPlot.m](./comparisonPlot.m).

## Caution
Due to the iterative process of updating concentration values, it has also been observed that complete soil saturation happens within certain time steps, i.e., there will be no changes in the figure after that, though the code will run until it reaches the final time.

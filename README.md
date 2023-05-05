# MATH-566-Final-Project

## Description 

This program solves the 2D-transient seepage problem using the finite difference scheme. The goal of this program is to simulate groundwater seepage flow to visualize the impact of the non-diagonal hydraulic conductivity tensor on seepage flow through anisotropic soil for various stratigraphic tilts with respect to the horizontal axis. The model simulates seepage flow for a domain (length 2.1 m and thickness 2 m) consisting of two soil types deposited in two layers stacked up together with the same length, which can also be used to simulate the flow for a single layer of soil considering the fact that hydraulic conductivity for both of the types is the same. In addition, the code is flexible to simulate other scenarios, e.g., soil layers stacked either horizontally or vertically, and the degree of tilt of the soil stratigraphic planes with respect to the horizontal axis. The boundary walls of the domain are impermeable, but to simulate the flow, several types of inlet(s) and outlet(s) combinations on boundary walls are considered as follows:

           type A: Horizontal Gradient(Midpoint left inlet, Midpoint right outlet & Endpoints right outlets)

           type B: Vertical Gradient(Midpoint bottom  inlet, Midpoint top outlet & Endpoints top outlet) 

           type C: One-to-one Combiantion (Bottom left inlet, Top right outlet & Top left inlet, Bottom right outlet)

           type D: Concrete Dam (Upstream lake inlet & Downstream lake outlet)

The governing equation, boundary conditions, mathematical formulation, and implementation of the numerical scheme are discussed in the file attached named [Mathematics_behind_the_Artifact.pdf](./Mathematics_behind_the_Artifact.pdf).

## Configure Instructions
To run the program, one needs to download the files in the code folder, then open the [main.m](./main.m) in MATLAB interface. The program requires the user to provide values for certain parameters to suit the details of the desired analysis, though for a few parameters a default value can be chosen. The considered geometry of the problem requires at least 10 horizontal and vertical nodes. This program is specific to SI units.

## Output
The program will generate three figures. The first figure will show the distribution of nodes according to chosen inlet and outlet combinations; the second figure will show the plot of flow velocity with equipotential lines; and the last figure will show the plot of contour lines for the pressure head.

## A little-known bug
It is possible for the code to become unstable when certain parameter adjustments are made. Code may become unstable and give erroneous results if the hydraulic conductivity of the soil is too small. Due to the iterative process of updating hydraulic conductivity values, it has also been observed that complete soil saturation happens within the first time step.

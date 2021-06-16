# Diffusion2D_Library
<img src = "../images/Cfield1.png" width="100" align ="right">
Solves the diffusion equation in 1 or 2 dimensions.

<img src="2D_Region.png" width="400" align = "right">
"Diffusion2D_Library" is a C# library for solving parabolic partial differential equations in 1 or 2 dimensions.  A representative region, &Omega;, in the cartesian plane, with boundary conditions specified on &delta;&Omega; and initial condition indicated, is shown to the right.  

The classic examples for parabolic partial differential equations are the heat and the diffusion equations.  The forms for both equations are shown in Eq 1 and the solutions are subject to an initial and boundary conditions.  The differences lie in the physical interpretation of the terms that make up the &nu; parameter and the types of boundary conditions imposed.  We will focus on solutions to the diffusion equation in this write-up.

***

|1| &delta;c/&delta;t - D(&delta;<sup>2</sup>c/&delta;x<sup>2</sup>) - f(x,t) = 0|
|-|--------------------------------------------------------------------------------|


where  x is a 1 or 2-dimensional vector and D represents the diffusivity of the diffusing species.


Subject to an initial condition:


|2| c(x,0) = I<sub>0</sub>(x)|
|-|--------------------------|


and boundary conditions:

in 1 dimension:

|3| c(0,t) = g<sub>0</sub>(t)|
|-|--------------------------|

|4| c(L,t) = g<sub>1</sub>(t)| 
|-|--------------------------|

 or 2 dimensions:

|5| c(0,y,t) = g<sub>0</sub>(t)|
|-|----------------------------|

|6| c(L,y,t) = g<sub>1</sub>(t)|
|-|----------------------------|

|7| c(x,0,t) = g<sub>2</sub>(t)|
|-|----------------------------|

|8| c(x,H,t) = g<sub>3</sub>(t)|
|-|----------------------------|
 
***
## Numerical Solvers
 The region, &Omega;, is discretized in space such that &Delta;x = &Delta;y and there are N<sub>x</sub> grid points in the x-direction and N<sub>y</sub> grid points in the y-direction.  Time is discretized into time-steps of &Delta;t duration.  Diffusion is calculated for N<sub>t</sub> time-steps.  The Fourier mesh number, &nu;, for the discretized space and time is given by equation 9.
 

|9| &nu; = D&Delta;t/&Delta;x<sup>2</sup>|
|-|--------------------------------------|

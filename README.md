# ROLLED_UP_QW_spectrum

This repository shows the codes which I used in my physics bachelor's thesis at Industrial University of Santander(UIS), Colombia.   The codes have written in python 3 using Jupyter Lab that you can use installing Anaconda application. 

My thesis was about to study the spectral properties of a shall donor in a semiconductor n-type which is rolled-up in a form to spiral. 

Structure:

## Models Semi_circles:

The spiral consists from semicircle wings of different radii Rs. The curvature radius initially keeps the value R0 for the polar coordinate 0 < phi < π, then at the point of connection of two adjacent semicircles the radius Rs is increased ∆R, while the center of curvature is displaced to the right in ∆R. As the angle phi increasing further reaches the value the curvature radius is increased again and the curvature center is displaced to the left. The curvature radius R_k and the arc length s_k at the beginning of a semicircle wing with the number k, are equal to and R_k = R0 + ∆R · k and, s_k = (k + 1) · π · (R0 + 0:5k · ∆r), respectively.

## Model Archimedian Spiral:

R(phi) = a+ b/2pi * phi 

We do a interpolation to get s in terms of phi. 




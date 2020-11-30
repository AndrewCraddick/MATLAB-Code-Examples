# MATLAB-Code-Examples
Examples of my code written to simulate electric and magnetic fields produced by various static charge/current distributions

All simulations are performed using Colomb's Law or Biot Savart's Law. This causes the coding to be far more difficult because I do not utilize symmetry inherent to the charge/current distribution. The benefit is that simulating distributions using these laws rather than those which rely on symmetry (Gauss's Law and Ampere's Law) is that I am develop an approach to simulate the distribution which works under circumstances where the distribution lacks symmetry and a computationally heavy brute force method must be applied.

For this reason, the code is frequently longer than would typically be needed to simulate relatively simple situations such as a perfect ring of charge where assumptions based on symmetry drastically simplify and reduce the complexity of the approach needed to accurately simulate the distribution.

Nevertheless, I was very satisfied with my results and was forced to think creatively to mathematically represent the ditributions and impliment those representations in code which was efficient and accurate without a long execution time.

In this repository there are simulations of:

- Electric field and electric potential due to a ring of charge
- Magnetic field of a physical dipole
- Isosurfaces of equipotential due to a uniformly charged sphere *
- Isosurfaces of equipotential due to a charged sphere in the presence of a separate charged particle *
- Magnetic field of a short and "long" solenoid (which is modeled as a helix using parametric equations) *
- Magnetic field of a helmholtz coil *

 (more advanced examples are denoted with an "*")
 
 Some examples contain more explanatory comments than others.

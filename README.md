![Header](https://github.com/alicialawjy/MATLAB-sway-structure/blob/main/Screenshots/Truss_Configuration.png)

# Modelling of a Sway Structure
The Finite Element (FE) Method is a very powerful mathematical technique used to solve complex systems of partial differential equations. When applied with respect to structural analysis of trusses and frames, the method:
1. Breaks the structure down into a limited and countable number of elements, 
2. Reformulates the governing differential euqations (equilibrium, constitutive and kinematic), and finally 
3. Assembles the elements to obtain the full structure again. 
This allows us to build up a matrix of simultaneous algebraic equations that can now be solved. 

Using this concept, these procedures are applied to the codes contained in the classes: 
- ‘TRUSS’ saved in **TRUSS.m**, and 
- ‘FRAMES’ saved in **FRAMES.m**, 

under the functions ‘assemble’ and ‘solve’ within these respective classes. 

Script **Radius.m** contains the function that calculates the maximum and minimum sway displacement (i.e. how much the structure has moved) for both Truss and Frame systems.

Script **Qn1_rcond.m** is used to determine the number of cross bracings (i.e diagonal elements) required to achieve a statically stable structure via the rcond value.
To give some background:
- 'rcond' is known as the reciprocal condition number. It is a scale-invariant measure of how close a given matrix is to the set of singular matrices. 
- Computationally, if the rcond is smaller than the machine precison, ε = 1 x 10-15, the matrix is considered badly conditioned and the structure is statically unstable.
- Aim: find the limit where the rcond > ε.

Script **CourseworkScript.m** is used to run through different possible design configurations for both Truss and Frame systems to identify the most optimal structure. It takes into consideration:
- Economy: material quantity and costs,
- Displacement: how much the system moves, and
- Factor of Safety: how safe the structure is.

To look at the final results and findings, see **Computational Report.pdf** for the full report.

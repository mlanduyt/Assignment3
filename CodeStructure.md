## Overview of FEA Code

![alt text](codestructure.png)

### Code Overview

Execution of code is performed via full_code_example_1.py.
This code executes the remaining files in the following structure.

1. pre_proces.py generates the mesh using either triangular or quadrilateral elements. This file returns node coordinates and element connections. 

2. solver.py executes the analysis part of the code. Solver.py applies the boundary conditions and solves using the following subscripts:
    2a. discretization.py determines shape functions and derivatives; determines gauss point locations
    2b. assemble_global.py loops over all elements, calls local_element.py
    2c. local_element.py determines internal forces and distributes loads. 

3. visualize.py performs the post-processing aspects of the code. Displaying meshes and visualizing displacement magnitudes. This script generates 2D_meshs and gifs. 
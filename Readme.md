# 3D Level Set Topology Optimization via Reaction diffusion method

This repository contains MATLAB code for performing 3D level set topology optimization using the reaction-diffusion method.

## Files

### 1. LBeamLS.m

This is the main script that sets up and runs the optimization for an L-beam example.

Key features:
- Defines material properties
- Sets up geometry and boundary conditions
- Specifies optimization parameters
- Calls the main optimization function `RDLSTO3D`
- Plots iteration history

### 2. RDLSTO3D.m

This file contains the core optimization loop and helper functions for the 3D level set topology optimization.

Key functions:
- `MeshGeneration`: Creates the 3D mesh
- `BoundaryConditionsImplementation`: Applies boundary conditions
- `BrickElementStiffnessCalculation`: Calculates element stiffness matrices
- `RD_T1T2Generator`: Generates matrices for the reaction-diffusion method

The main optimization loop updates the level set function using the reaction-diffusion method and calculates compliance and volume fractions.

### 3. PlotI3DterationHistory.m

This function visualizes the optimization results.

Features:
- Plots the 3D structure at each iteration
- Allows customization of the view and axis order
- Uses isosurfaces to represent the optimized structure

## Usage

1. Set up your problem parameters in `LBeamLS.m`
2. Run `LBeamLS.m` to perform the optimization
3. The results will be visualized using `PlotI3DterationHistory`

## Dependencies

This code requires MATLAB with the following toolboxes:
- Partial Differential Equation Toolbox
- Image Processing Toolbox

## Note

This code is designed for research and educational purposes. It may require further optimization for large-scale industrial applications.

## Contributing

Contributions to improve the code or extend its capabilities are welcome. Please submit pull requests or open issues for any bugs or feature requests.

## Results


https://github.com/user-attachments/assets/a6c8cb67-6dae-4b7a-8f34-581257527c2a


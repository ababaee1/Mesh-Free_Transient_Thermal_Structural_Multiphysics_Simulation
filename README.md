# Mesh-Free Transient Thermal/Structural Multi-Physics Simulation

This MATLAB repository contains a **Mesh Free transient (quasi-static) thermal/structural simulation** model designed for multi-physics analysis. The model simulates the interaction of thermal and structural fields under transient thermal loading conditions and outputs deformation, stress variants, and temperature distribution over time and space.

## Key Features

- **Nonlinear Analysis**: Uses the Newton-Raphson algorithm to handle nonlinear equations, critical for accurate structural simulation under thermal loads.
- **Generalized Differential Quadrature (GDQ)**: Employed for spatial discretization, allowing mesh free efficient and high-accuracy solution over the domain.
- **Boundary Conditions**: The model includes varied boundary conditions, such as insulation and prescribed temperatures, to simulate realistic thermal interactions.
- **Temporal and Spatial Resolution**: Outputs simulation results over both time and spatial domains, providing detailed insight into deformation and stress distribution.

## Outputs

1. **Deformation**: Spatial and temporal deformation profiles of the structure.
2. **Stress Variants**: Calculates stress distribution resulting from thermal gradients and structural constraints.
3. **Temperature Profile**: Detailed temperature distribution over time.

## Parameters

The model is set up with the following configurable parameters:

- **nr**: Number of nodes in the radial direction (default: 15).
- **nz**: Number of nodes in the z-direction (default: 15).
- **dt**: Time step for the heat conduction equation (default: 0.01).
- **Material Dimensions**: Inner and outer radius, and plate thickness (e.g., `a=100e-3`, `b=1e-8`, `h=5e-3`).
- **Boundary Conditions**: Options include:
  - **`insul`** for zero heat flux
  - **`convection`** for heat convection with air
  - **`prescribed`** for a specific temperature setting
- **Partial Load Ratios**: Configurable load distribution on the top and bottom surfaces.

## Code Usage

1. **Initialize Parameters**: Update the desired values in the main script to set up the grid and boundary conditions.
2. **Run Simulation**: Execute the main script. The simulation will compute the deformation, temperature, and stress fields.
3. **Plot Results**: Optional plotting is available at the end of the script to visualize the temperature and deformation profiles.

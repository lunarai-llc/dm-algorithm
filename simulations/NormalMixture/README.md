# Two-Component Normal Mixture Simulation

This folder contains simulation code for a **two-component normal mixture model**
used in the numerical experiments of the accompanying paper.

## Model specification

The true data-generating model is a two-component normal mixture with:

- Mixing proportion: π₁ = 0.3
- Component 1: μ₁ = 10, σ₁ = 1
- Component 2: μ₂ = 0, σ₂ = 1
- Sample size: n = 200

## Simulation setup

- The number of mixture components **K is unknown** and is estimated using the
  **Divergence Information Criterion (DIC)**.
- **Data splitting is not used** in this simulation.
- Density estimation is performed using an **Epanechnikov kernel**.

## Main script

- `simulation_normal_mixture.m`  
  Runs the Monte Carlo simulation for the setting described above.

## Output and paper reference

The results produced by this simulation correspond to the table titled:

> *“Two-Component Normal (Epanechnikov kernel; 100% chose K = 2)”*

in the accompanying paper.

## Computational notes

Default settings are configured for reproducibility and moderate runtime.
Full-scale Monte Carlo settings used in the paper are documented in the script
comments.

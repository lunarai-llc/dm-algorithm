# Poisson Mixture Simulations

This folder contains simulation and figure-generation code for
**two-component Poisson mixture models** used in the numerical experiments
of the accompanying paper.

The true data-generating model throughout is a two-component Poisson mixture:
\[
0.4\,\text{Poi}(0.5) + 0.6\,\text{Poi}(10).
\]

---

## Model specification

- Mixing proportion: π₁ = 0.4
- Component 1: Poisson(λ₁ = 0.5)
- Component 2: Poisson(λ₂ = 10)
- Sample size: n varies depending on the experiment
- Number of mixture components: K = 2 (known)

---

## Prerequisites

Before running any simulation in this folder, ensure that all helper functions
are available by adding the `src/` directory to the MATLAB path:

```matlab
addpath(genpath("src"));

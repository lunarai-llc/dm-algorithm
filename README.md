# Divergence-Minimization (DM) Algorithm for Latent-Structure Models

This repository provides a **reference MATLAB implementation** of the
**Divergence-Minimization (DM) algorithm** for latent-structure models, as proposed in:

> **Lei Li and Anand N. Vidyashankar (2026).**  
> *Divergence-Minimization for Latent-Structure Models:  
> Monotone Operators, Contraction Guarantees, and Robust Inference.*  

The DM framework generalizes EM by replacing likelihood optimization with
divergence minimization, yielding **monotone descent**, **local contraction**, and
**robust inference** under bounded residual-adjustment divergences (e.g.,
Hellinger, negative exponential).

---

## Scope of this repository

This codebase is intended for **research reproducibility** and accompanies the paper.
It includes:

- A reference implementation of the DM algorithm for finite mixture models
- Instantiations using Hellinger and negative-exponential divergences
- Simulation scripts to reproduce representative numerical results

The implementation prioritizes **clarity and correctness** over computational
optimization.

---

## Directory structure

```text
dm-algorithm/
├── src/            # Core DM algorithm and divergence utilities
├── simulations/    # Simulation drivers reproducing paper results
├── figures/        # Figure-generation scripts
├── output/         # Generated results (ignored by git)
├── README.md
└── LICENSE

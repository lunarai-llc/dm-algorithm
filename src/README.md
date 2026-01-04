# Core Functions (`src/`)

This folder contains **all necessary MATLAB functions** required to reproduce
the simulation results and data analyses in the accompanying paper.

---

## Nonparametric density estimation (kernels)

- `Bino_kernel.m`  
  Function to calculate the nonparametric estimate using the Binomial kernel.

- `Poisson_kernel.m`  
  Function to calculate the nonparametric estimate using the Poisson kernel.

- `NB_kernel.m`  
  Function to calculate the nonparametric estimate using the Negative Binomial kernel.

- `Example1_kernel.m`  
  Function to calculate the triangle-kernel-based estimate.

- `tri_zero_mass_bandwidth.m`  
  Function for discrete triangular kernel estimation with zero-mass control for estimating nonparametric \( g_n \).

---

## Data generation functions

### Poisson–Lognormal (PL) mixture models
- `Generate_PL1.m` – generate 1-component Poisson–Lognormal model data  
- `Generate_PL2.m` – generate 2-component Poisson–Lognormal model data  
- `Generate_PL3.m` – generate 3-component Poisson–Lognormal model data  
- `Generate_PL4.m` – generate 4-component Poisson–Lognormal model data  

### Poisson mixture models
- `Generate_Poi2.m` – generate 2-component Poisson mixture model data  

### Poisson–Gamma (PG) mixture models
- `Generate_Y_iid_PG2.m` – generate 2-component PG model data  
- `Generate_Y_iid_PG3.m` – generate 3-component PG model data  
- `Generate_Y_iid_PG4.m` – generate 4-component PG model data  
- `Generate_Y_iid_PG5.m` – generate 5-component PG model data  

---

## EM algorithm

- `EM_iid_Poi2.m`  
  EM algorithm for a two-component Poisson mixture model.

- `EM_iid_PG2.m`  
  EM algorithm for a two-component Poisson–Gamma mixture model.

- `EM_iid_PG2_HD.m`  
  HMIX (Hellinger-distance-based) algorithm for a two-component PG model.

- `EM_iid_PG2_NED.m`  
  VNEDMIX algorithm for a two-component PG model.

---

## HMIX algorithms

### Two-component Poisson mixtures
- `HMIX_two_Poisson_Bino_kernel.m` – Binomial kernel  
- `HMIX_two_Poisson_Empirical_kernel.m` – empirical kernel  
- `HMIX_two_Poisson_Example1_kernel.m` – triangular kernel  
- `HMIX_two_Poisson_NB_kernel.m` – Negative Binomial kernel  
- `HMIX_two_Poisson_Poisson_kernel.m` – Poisson kernel  

### Two-component Poisson–Gamma mixtures
- `HMIX_PG2_Bino_kernel.m` – Binomial kernel  
- `HMIX_PG2_Ex1_kernel.m` – triangular kernel (c = 1)  
- `HMIX_PG2_NB_kernel.m` – Negative Binomial kernel  
- `HMIX_PG2_Poiss_kernel.m` – Poisson kernel  

---

## NEDMIX / VNEDMIX algorithms

### Two-component Poisson mixtures
- `NEDMIX_two_Poisson_Bino_kernel.m` – Binomial kernel  
- `NEDMIX_two_Poisson_Empirical_kernel.m` – empirical kernel  
- `NEDMIX_two_Poisson_Example1_kernel.m` – triangular kernel  
- `NEDMIX_two_Poisson_NB_kernel.m` – Negative Binomial kernel  
- `NEDMIX_two_Poisson_Poisson_kernel.m` – Poisson kernel  

### Two-component Poisson–Gamma mixtures
- `NEDMIX_PG2_Bino_kernel.m` – Binomial kernel  
- `NEDMIX_PG2_Ex1_kernel.m` – triangular kernel  
- `NEDMIX_PG2_NB_kernel.m` – Negative Binomial kernel  
- `NEDMIX_PG2_Poiss_kernel.m` – Poisson kernel  

---

## Initialization routines

- `INI_Poi2.m`  
  K-means-based initialization for two-component Poisson mixture models.

- `INI_PG2.m` – initialization for 2-component PG model  
- `INI_PG3.m` – initialization for 3-component PG model  
- `INI_PG4.m` – initialization for 4-component PG model  
- `INI_PG5.m` – initialization for 5-component PG model  

- `IID_PG_solve_alpha_beta.m`  
  Moment-based initial value estimation for PG model parameters.

---

## Integrated Squared Error (ISE)

### Poisson mixture models
- `ISE_two_Poisson_Bino_kernel.m` – Binomial kernel  
- `ISE_two_Poisson_Empirical.m` – empirical kernel  
- `ISE_two_Poisson_Example1.m` – triangular kernel  
- `ISE_two_Poisson_NB_kernel.m` – Negative Binomial kernel  
- `ISE_two_Poisson_Poisson_kernel.m` – Poisson kernel  

### Poisson–Gamma mixture models
- `ISE_two_PG_Bino_kernel.m` – Binomial kernel  
- `ISE_two_PG_Empirical.m` – empirical kernel  
- `ISE_two_PG_Example1.m` – triangular kernel (c = 1)  
- `ISE_two_PG_NB_kernel.m` – Negative Binomial kernel  
- `ISE_two_PG_Poisson_kernel.m` – Poisson kernel  

---

## Model selection and data splitting

- `selectK_DIC.m`  
  Function to select the number of mixture components \( K \) using the DIC criterion.

- `selectK_normmix.m`  
  Function to estimate the number of components \( K \) and parameter estimates
  for normal mixture models using EM, HMIX, and VNEDMIX with Epanechnikov kernel
  estimation of \( g_n \).

- `mix_split_select_estimate.m`  
  Data-splitting procedure for estimating unknown \( K \) using EM, HMIX, and
  VNEDMIX for Poisson, PG, and PL mixture models.

- `mix_mc_identifyK.m`  
  Complete data-splitting algorithm for identifying unknown \( K \) with variance
  estimation.

- `mix_mc_identifyK_perMethodFast.m`  
  Faster variant of the data-splitting algorithm without variance computation.

- `summarize_pg_results_data_splitting.m`  
  Summarize the table titled *“Parameter Estimation with Unknown K with Data
  Splitting Method”* for PG mixture models.

---

## Unified interface

- `Mixture_Model_No_Regress.m`  
  A unified function implementing **EM, HMIX, and NEDMIX** algorithms for finite
  mixture models with an arbitrary number of components.

  An example usage of this function is provided in:

---

## Demonstration and plotting utilities

- `demo_normmix_epa_sim.m`  
Simulation for a two-component normal mixture with unknown \( K \)
(π₁ = 0.3, μ₁ = 10, σ₁ = 1, μ₂ = 0, σ₂ = 1, n = 200).

- `demo_pg_outlier_sim.m`  
PG mixture simulation with contamination levels ranging from 0% to 30%
(π₁ = 0.3, α₁ = 10, β₁ = 1, α₂ = 1, β₂ = 2).

- `plot_2Pois_figure.m`  
Plot of average parameter estimates using EM, HMIX, and VNEDMIX across sample
sizes. Corresponds to the table titled
*“Estimates Using Different Kernels for 0.4 Poi(0.5) + 0.6 Poi(10)”*.

---

## Image segmentation utilities

- `em_poissmix_bic.m` – EM-based Poisson mixture image segmentation  
- `em_bic_plot.m` – EM-based image reconstruction and plotting  
- `hmix_hic_plot.m` – HMIX-based image segmentation  
- `vnedmix_dic_plot.m` – VNEDMIX-based image segmentation  

---

## Notes

- This folder contains **functions only**; no simulation outputs or data files
are stored here.
- All simulation and figure scripts assume that this folder has been added to
the MATLAB path:
```matlab
addpath(genpath("src"));

---


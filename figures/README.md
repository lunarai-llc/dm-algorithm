# Figure Generation (`figures/`)

This folder contains datasets and scripts used to generate figures for the
**image segmentation analysis** based on Poisson mixture models.

The scripts apply **EM, HMIX, and VNEDMIX algorithms** to both clean and
contaminated image data and are used to produce figures reported in the paper.

---

## Datasets

- `Image_a.mat`  
  Phantom image dataset in MATLAB format.

- `Image_c.mat`  
  Lena image dataset in MATLAB format.

- `image_c_contaminated_pois250_prob_3.mat`  
  Lena image dataset with 30% contamination outliers.

---

## Figure-generation scripts

- `image_a_phantom_image_analysis.m`  
  Generate figures by analyzing the phantom image using EM, HMIX, and VNEDMIX
  algorithms under a Poisson mixture model (no contamination).

- `example_image_c_lena_analysis.m`  
  Generate figures by analyzing the Lena image using EM, HMIX, and VNEDMIX
  algorithms under a Poisson mixture model (no contamination).

- `example_image_c_lena_analysis_contaminated30.m`  
  Generate figures by analyzing the Lena image using EM, HMIX, and VNEDMIX
  algorithms under a Poisson mixture model with 30% contamination.

---

## Notes

- All scripts assume that the `src/` folder has been added to the MATLAB path:
  ```matlab
  addpath(genpath("src"));

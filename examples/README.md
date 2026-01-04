# Example Data Analysis (`examples/`)

This folder contains example datasets and analysis scripts used to demonstrate
the application of the proposed methods to **image segmentation problems**
using Poisson mixture models.

The examples illustrate the use of **EM, HMIX, and VNEDMIX algorithms** on both
clean and contaminated image data.

---

## Example datasets

- `Image_a.mat`  
  Phantom image dataset in MATLAB format.

- `Image_c.mat`  
  Lena image dataset in MATLAB format.

- `image_c_contaminated_pois250_prob_3.mat`  
  Lena image dataset with 30% contamination, where outliers are generated from
  a Poisson distribution with mean 250.

---

## Example analysis scripts

- `image_a_phantom_image_analysis.m`  
  Analysis of the phantom image using EM, HMIX, and VNEDMIX algorithms under a
  Poisson mixture model (no contamination).

- `example_image_c_lena_analysis.m`  
  Analysis of the Lena image using EM, HMIX, and VNEDMIX algorithms under a
  Poisson mixture model (no contamination).

- `example_image_c_lena_analysis_contaminated30.m`  
  Analysis of the Lena image using EM, HMIX, and VNEDMIX algorithms under a
  Poisson mixture model with 30% contamination.

---

## Notes

- All example scripts assume that the `src/` folder has been added to the
  MATLAB path:
  ```matlab
  addpath(genpath("src"));

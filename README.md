# Lamian Helper Package

This package provides a convenient interface to use **Lamian**, a tool designed to compare expression changes over pseudotime across different conditions. For detailed information on Lamian, please refer to Professor Wenpin Zhou's original repository:  
https://github.com/Winnie09/Lamian/tree/master

---

## Features

Starting from a Seurat object containing **sample** and **gene cluster (meta.data)** information, this package enables the following workflow:  

1. **Harmony integration & normalization**  
2. **Subsetting & renaming samples or metadata**  
3. **Slingshot trajectory inference analysis**  
4. **Filtering of trajectory and features**  
5. **Lamian parameter configuration**  
6. **Visualization options:**  
   - Heatmap  
   - Population plot  
   - FeaturePlot split by sample  

---

## Job Submission

- Tasks are submitted via `qsub` with **28 CPU cores** allocated.  
- Integrated **email notification** functionality for job status updates.

---

## Important Fixes & Modifications

We have incorporated fixes and improvements in the Lamian codebase to ensure robustness and usability:

1. **Avoiding no-solution problems** during trajectory fitting  
2. **Fixing matrix non-invertibility issues**  
3. **Adjusting some default parameters and visualization settings** for better performance and clarity  

These fixes are essential to prevent errors during analysis.

---

## Installation

The corrected version of Lamian and all necessary code reside in the `Lamian` folder of this repository. To install, use `devtools`:  

```r
devtools::install_local("path_to_this_repo/Lamian")
```

---

## Example Usage

The typical workflow includes loading your Seurat object, running normalization, harmony, slingshot, and then using Lamian for comparative pseudotime analysis with the enhanced visualization options provided.

---

## Contact

For any questions or issues, please feel free to open an issue in this repository or contact the maintainer.

---

**Note:** Make sure to use this modified version of Lamian to avoid runtime errors and take advantage of the enhanced stability and visualization features.
